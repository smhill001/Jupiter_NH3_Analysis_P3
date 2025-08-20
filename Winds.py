import numpy as np
import csv

try:
    from scipy.signal import correlate2d
except ImportError:
    raise ImportError("scipy is required for this script. Please install or update scipy via pip or conda.")

JUPITER_RADIUS_M = 71492e3

# ----------------------- Mask (optional) -----------------------

def load_mask(mask_path, shape):
    from imageio import imread
    if mask_path is None:
        return None
    mask = imread(mask_path, as_gray=True)
    if mask is None:
        raise ValueError(f"Could not read mask at {mask_path}")
    mask = (mask > 0).astype(np.float32)
    if mask.shape != shape:
        from skimage.transform import resize
        mask = resize(mask, shape, order=0, preserve_range=True, anti_aliasing=False).astype(np.float32)
    return mask

# ------------------ Correlation-based shift --------------------

def _parabolic_subpixel_1d(fm, f0, fp):
    denom = max((fm - 2*f0 + fp), 1e-12)
    return 0.5 * (fm - fp) / denom

def correlate_shift(img1, img2, use_window=True, normalize=True, boundary='symm'):
    """Estimate (dx, dy) that moves img1 -> img2 using normalized, windowed xcorr.

    Returns subpixel shifts using simple quadratic peak interpolation.
    """
    a = np.ascontiguousarray(img1, dtype=np.float32)
    b = np.ascontiguousarray(img2, dtype=np.float32)
    if a.ndim != 2 or b.ndim != 2:
        raise ValueError("correlate_shift requires 2D grayscale input")

    # reject flat tiles early
    if a.std() < 1e-6 or b.std() < 1e-6:
        return 0.0, 0.0, 0.0, 0.0

    if use_window:
        wy = np.hanning(a.shape[0]).astype(np.float32)
        wx = np.hanning(a.shape[1]).astype(np.float32)
        w = np.outer(wy, wx)
        a = a * w
        b = b * w

    if normalize:
        a = (a - a.mean()) / (a.std() + 1e-6)
        b = (b - b.mean()) / (b.std() + 1e-6)

    # Full correlation avoids center/offset ambiguity; boundary='symm' reduces edge bias
    corr = correlate2d(b, a, mode='full', boundary=boundary)
    y0, x0 = np.unravel_index(np.argmax(corr), corr.shape)
    dy = y0 - (a.shape[0] - 1)
    dx = x0 - (a.shape[1] - 1)

    # Peak SNR proxy: peak vs. 99th percentile of remainder
    peak = corr[y0, x0]
    corr_flat = corr.ravel()
    corr_flat[y0 * corr.shape[1] + x0] = np.min(corr_flat)
    ref = np.percentile(corr_flat, 99.0)
    snr = float(peak / (ref + 1e-6))

    # Sub-pixel refinement (parabolic fit) if not on border
    if 0 < y0 < corr.shape[0]-1 and 0 < x0 < corr.shape[1]-1:
        dx_sub = _parabolic_subpixel_1d(corr[y0, x0-1], corr[y0, x0], corr[y0, x0+1])
        dy_sub = _parabolic_subpixel_1d(corr[y0-1, x0], corr[y0, x0], corr[y0+1, x0])
        dx += dx_sub
        dy += dy_sub

    return float(dx), float(dy), float(peak), snr

# ---------------- Tile correlation (sparse) --------------------

def tile_correlation(img1, img2, tile, stride, min_snr=1.1):
    h, w = img1.shape
    vectors = []
    for y in range(0, h - tile + 1, stride):
        for x in range(0, w - tile + 1, stride):
            roi1 = np.ascontiguousarray(img1[y:y + tile, x:x + tile], dtype=np.float32)
            roi2 = np.ascontiguousarray(img2[y:y + tile, x:x + tile], dtype=np.float32)
            dx, dy, peak, snr = correlate_shift(roi1, roi2)
            if snr < min_snr:
                continue  # reject flat/ambiguous tiles
            vectors.append((x + tile / 2.0, y + tile / 2.0, dx, dy))
    return vectors

# ----------------- Farnebäck (sampled grid) --------------------

def _to_u8(img):
    # robust contrast stretch to uint8
    lo, hi = np.percentile(img, (1, 99))
    if hi <= lo:
        hi = lo + 1.0
    out = np.clip((img - lo) * (255.0 / (hi - lo)), 0, 255).astype(np.uint8)
    return out

FARNEBACK_PRESETS = {
    # Smooth, large-scale jets over full 1800x3600 maps (good first pass)
    "jets_global": dict(pyr_scale=0.5, levels=5, winsize=25, iterations=5,
                         poly_n=7, poly_sigma=1.5, flags=0),
    # Mixed scales (belts/eddies) with moderate detail
    "eddies_medium": dict(pyr_scale=0.5, levels=4, winsize=17, iterations=4,
                           poly_n=5, poly_sigma=1.2, flags=0),
    # Fine features (spots, small vortices); smaller displacements
    "fine_features": dict(pyr_scale=0.5, levels=3, winsize=11, iterations=3,
                           poly_n=5, poly_sigma=1.1, flags=0),
}

def farneback_flow(img1, img2, step=12, preset=None, **kwargs):
    import cv2
    # Preprocess to stable uint8 with robust contrast
    a = _to_u8(img1)
    b = _to_u8(img2)

    # Optional light blur to suppress pixel noise without killing bands
    pre_blur_sigma = kwargs.pop("pre_blur_sigma", 0.0)
    if pre_blur_sigma and pre_blur_sigma > 0:
        k = max(3, int(pre_blur_sigma * 3) * 2 + 1)
        a = cv2.GaussianBlur(a, (k, k), pre_blur_sigma)
        b = cv2.GaussianBlur(b, (k, k), pre_blur_sigma)

    # Equalize histogram (optional)
    if kwargs.pop("equalize", False):
        a = cv2.equalizeHist(a)
        b = cv2.equalizeHist(b)

    params = FARNEBACK_PRESETS.get(preset, None)
    if params is None:
        # fallback to kwargs or sensible defaults
        params = dict(pyr_scale=0.5, levels=4, winsize=17, iterations=4,
                      poly_n=5, poly_sigma=1.2, flags=0)
    # allow caller to override any preset value via kwargs
    params.update(kwargs)

    flow = cv2.calcOpticalFlowFarneback(a, b, None, **params)
    h, w = a.shape
    vectors = []
    for y in range(0, h, step):
        for x in range(0, w, step):
            dx, dy = flow[y, x]
            vectors.append((float(x), float(y), float(dx), float(dy)))
    return vectors

# ----------------- Geo conversion helpers ----------------------


def deg_to_meters_per_px(deg_per_px, ref_latitude_deg):
    lat_rad = np.deg2rad(ref_latitude_deg)
    meters_per_deg_lon = (2 * np.pi * JUPITER_RADIUS_M * np.cos(lat_rad)) / 360.0
    return deg_per_px * meters_per_deg_lon

# -------------------------- Main --------------------------------

def estimate_winds(
    img1_path,
    img2_path,
    dt_seconds,
    meters_per_pixel=None,
    deg_per_px=None,
    ref_latitude_deg=0.0,
    #!!Does this set a single value for conversion of deg long to meters?
    method="phasecorr",
    tile=65,              # odd size avoids center ambiguity
    stride=48,
    mask_path=None,
    out_csv_path=None,
    out_quiver_path=None,
    farneback_params=None,
    farneback_step=8,
    min_snr=1.1
):
    from imageio import imread
    img1 = np.ascontiguousarray(imread(img1_path, as_gray=True).astype(np.float32))
    img2 = np.ascontiguousarray(imread(img2_path, as_gray=True).astype(np.float32))

    if meters_per_pixel is None and deg_per_px is not None:
        meters_per_pixel = deg_to_meters_per_px(deg_per_px, ref_latitude_deg)
    elif meters_per_pixel is None:
        raise ValueError("Must provide either meters_per_pixel or deg_per_px.")

    mask = None
    if mask_path is not None:
        mask = load_mask(mask_path, img1.shape)

    if method == "phasecorr":
        vectors = tile_correlation(img1, img2, tile, stride, min_snr=min_snr)
    elif method == "farneback":
        fb_params = farneback_params or dict(
            pyr_scale=0.5, levels=3, winsize=15, iterations=3,
            poly_n=5, poly_sigma=1.2, flags=0
        )
        vectors = farneback_flow(img1, img2, step=farneback_step, **fb_params)
    else:
        raise ValueError("Unknown method")

    results = []
    for x, y, dx, dy in vectors:
        speed_m_s = (np.hypot(dx, dy) * meters_per_pixel) / dt_seconds
        direction_deg = (np.degrees(np.arctan2(-dx, -dy)) + 360) % 360
        results.append((x, y, dx, dy, speed_m_s, direction_deg))

    if out_csv_path:
        with open(out_csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["x", "y", "dx_px", "dy_px", "speed_m_s", "dir_from_deg"])
            writer.writerows(results)

    if out_quiver_path:
        import matplotlib.pyplot as plt
        xs = np.array([r[0] for r in results])
        ys = np.array([r[1] for r in results])
        dxs = np.array([r[2] for r in results])-0.
        dys = np.array([r[3] for r in results])
        # thin quiver automatically if too dense
        max_arrows = 10000
        if xs.size > max_arrows:
            sel = np.linspace(0, xs.size-1, max_arrows).astype(int)
            xs, ys, dxs, dys = xs[sel], ys[sel], dxs[sel], dys[sel]
        import matplotlib.pyplot as plt
        fig,axs=plt.subplots(1,figsize=(10.0,5.0), dpi=150, facecolor="black")
        axs.imshow(img1, cmap="gray")
        axs.quiver(xs, ys, dxs, dys, color="red", angles="xy", scale_units="xy", scale=0.3)

        plt.figure(figsize=(10, 5))
        plt.imshow(img1, cmap="gray")
        plt.quiver(xs, ys, dxs, dys, color="red", angles="xy", scale_units="xy", scale=0.3)
        plt.savefig(out_quiver_path, dpi=150)
        plt.close()

    return results

# -------------------- 1D Zonal Winds (per latitude) --------------------

def _latitudes_from_height(h):
    """Infer planetographic latitude centers for a global cylindrical map.
    Assumes rows cover +90° (north) to -90° (south) uniformly.
    Returns array of length h with latitude at pixel-row centers.
    """
    return 90.0 - (np.arange(h) + 0.5) * (180.0 / h)

def _meters_per_px_lon_at_lat(lat_deg, meters_per_pixel, deg_per_px, ref_latitude_deg):
    """Compute longitudinal meters-per-pixel at a given latitude.
    Priority: use deg_per_px if provided (physically correct). Otherwise, if
    meters_per_pixel is provided at ref_latitude, scale by cos(lat)/cos(ref_lat).
    """
    lat_rad = np.deg2rad(lat_deg)
    if deg_per_px is not None:
        meters_per_deg_lon = (2 * np.pi * JUPITER_RADIUS_M * np.cos(lat_rad)) / 360.0
        return deg_per_px * meters_per_deg_lon
    # fallback: scale provided meters_per_pixel by cos(latitude)
    ref_rad = np.deg2rad(ref_latitude_deg)
    scale = np.cos(lat_rad) / max(np.cos(ref_rad), 1e-6)
    return meters_per_pixel * scale

def _xcorr1d_circular(s1, s2):
    """Circular cross-correlation (s1 -> s2) via FFT. Returns subpixel shift in px.
    Positive shift means s1 must move +x (east/right) to match s2.
    """
    n = s1.size
    a = np.asarray(s1, dtype=np.float32)
    b = np.asarray(s2, dtype=np.float32)

    # normalize to zero-mean, unit-variance to suppress DC bias
    a = (a - a.mean()) / (a.std() + 1e-6)
    b = (b - b.mean()) / (b.std() + 1e-6)

    fa = np.fft.rfft(a)
    fb = np.fft.rfft(b)
    r = np.fft.irfft(fb * np.conj(fa), n=n)

    k = int(np.argmax(r))
    # quad subpixel around peak (circular neighbors)
    km = (k - 1) % n
    kp = (k + 1) % n
    denom = max((r[km] - 2*r[k] + r[kp]), 1e-12)
    delta = 0.5 * (r[km] - r[kp]) / denom  # in [-0.5, 0.5] typically
    shift = k + float(np.clip(delta, -0.5, 0.5))

    # unwrap to signed shift in [-n/2, n/2)
    if shift > n/2:
        shift -= n

    # simple SNR proxy
    peak = r[k]
    ref = np.partition(r, -10)[-10] if n >= 10 else np.max(r)
    snr = float(peak / (ref + 1e-6))

    return float(shift), float(peak), snr


def estimate_zonal_winds(
    img1_path,
    img2_path,
    dt_seconds,
    meters_per_pixel=None,
    deg_per_px=None,
    ref_latitude_deg=0.0,
    out_csv_path=None
):
    """
    Compute a single zonal (east-west) velocity at each latitude by
    performing 1D circular cross-correlation of each latitude row across all longitudes.

    Inputs match Cloud Advection Winds primary inputs. If `deg_per_px` is given, it
    is interpreted as longitudinal degrees-per-pixel; meters-per-pixel in longitude
    is then computed per-latitude using cos(lat). If only `meters_per_pixel` is given,
    it is taken at `ref_latitude_deg` and scaled ∝ cos(lat).

    Returns: lat_deg (N,), u_m_s (N,) with eastward positive.
    
    NOTE: THIS CODE ONLY WORKS FOR FULL LATITUDE RANGE = 90N - 90S
    """
    from imageio import imread
    import pylab as pl

    im1 = np.ascontiguousarray(imread(img1_path, as_gray=True).astype(np.float32))
    im2 = np.ascontiguousarray(imread(img2_path, as_gray=True).astype(np.float32))
    if im1.shape != im2.shape:
        raise ValueError("Images must have identical shape for zonal winds.")

    h, w = im1.shape
    lat_deg = _latitudes_from_height(h)

    u_m_s = np.full(h, np.nan, dtype=np.float32)
    snr_arr = np.full(h, np.nan, dtype=np.float32)

    for y in range(h):
        s1 = im1[y, :]
        s2 = im2[y, :]
        # skip near-blank rows
        if s1.std() < 1e-6 or s2.std() < 1e-6:
            continue
        shift_px, peak, snr = _xcorr1d_circular(s1, s2)
        mpp = _meters_per_px_lon_at_lat(lat_deg[y], meters_per_pixel, deg_per_px, ref_latitude_deg)
        u_m_s[y] = (shift_px * mpp) / dt_seconds  # eastward positive
        snr_arr[y] = snr

    # Optional: light smoothing in latitude to stabilize profiles
    # (Savitzky–Golay-like 1-2-1 kernel applied once)
    kernel = np.array([1, 2, 1], dtype=np.float32)
    kernel = kernel / kernel.sum()
    u_sm = np.copy(u_m_s)
    u_sm[1:-1] = (u_m_s[:-2] * kernel[0] + u_m_s[1:-1] * kernel[1] + u_m_s[2:] * kernel[2])

    # replace gross outliers (>5× MAD) with smoothed values
    med = np.nanmedian(u_sm)
    mad = np.nanmedian(np.abs(u_sm - med)) + 1e-6
    bad = np.abs(u_sm - med) > (5.0 * 1.4826 * mad)
    u_sm[bad] = med

    if out_csv_path:
        with open(out_csv_path, "w", newline="") as f:
            wtr = csv.writer(f)
            wtr.writerow(["latitude_deg", "u_east_m_s", "snr"])
            for lat, u, s in zip(lat_deg, u_sm, snr_arr):
                wtr.writerow([lat, u, s])

    fig,axs=pl.subplots(1,figsize=(8,6), dpi=150, facecolor="white")
    axs.plot(lat_deg,u_sm)
    return lat_deg, u_sm
