import numpy as np
import csv

JUPITER_RADIUS_M = 71492e3

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
    im1,
    im2,
    dt_seconds,
    meters_per_pixel=None,
    deg_per_px=None,
    ref_latitude_deg=0.0,
    out_csv_path=None):
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
    import pylab as pl
    import numpy as np

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

    return lat_deg, u_sm
