# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 12:57:12 2025

@author: smhil
"""
import numpy as np
from skimage.feature import peak_local_max
from astropy.time import Time
from copy import deepcopy
import csv


def convert_time(row_inds, col_inds, mean_time_array):
    fits_times = []
    jd_times = []
    print("%%%%%%%%%%row_inds, col_inds=",row_inds, col_inds)
    #print("mean_time_array[row_inds[0], col_inds[0]]=",mean_time_array[row_inds[0], col_inds[0]])
    for row, col in zip(row_inds, col_inds):
        try:
            #print("mean_time_array[row, col]=",mean_time_array[row, col],row,col)
            #temptime = Time(mean_time_array[row, col], format='jd')
            temptime = Time(mean_time_array[row, col], format='jd')
            fits_times.append(temptime.fits)
            tempJD = deepcopy(temptime)
            tempJD.format = 'jd'
            jd_times.append(tempJD.value)
        except:
            fits_times.append("N/A")
            jd_times.append(0.0)
    return fits_times, jd_times


def find_blob(image_to_segment, intensity_image,threshold_abs=None, mode='max'):
    import numpy as np
    import pylab as pl
    from skimage.measure import label, regionprops
    from skimage.io import imread, imshow

    if threshold_abs is None:
        threshold_abs = np.nanmean(image_to_segment)

    if mode == 'min':
        blob_mask =  image_to_segment < threshold_abs
    else:
        blob_mask = image_to_segment > threshold_abs
        
    labeled_image = label(blob_mask)
    props_data=regionprops(labeled_image,intensity_image=image_to_segment)
    props_intensity=regionprops(labeled_image,intensity_image=intensity_image)
    #for region in regionprops(labeled_image,intensity_image=intensity_image):
    #    print(f"Blob Label: {region.label}, Area: {region.area}, Centroid: {region.centroid}, \
    #          Centroid_Weighted: {region.centroid_weighted}")        
        
    return blob_mask,labeled_image,props_data,props_intensity


def process_blob(image_to_segment, intensity_image, lats, lon_lims, threshold_abs=None, mode='max'):
    """
    Segments blobs, converts coordinates to lat/lon, merges regionprops from both images,
    sorts by longitude, and relabels regions accordingly.

    Parameters:
        image_to_segment (ndarray): Image used for segmentation.
        intensity_image (ndarray): Image used for intensity-based measurement.
        lats (ndarray): Latitude array (1D or 2D) from the image grid.
        lon_lims (tuple): (min_lon, max_lon) representing the longitude extent.
        threshold_abs (float): Threshold for segmentation.
        mode (str): 'max' or 'min' for segmentation.

    Returns:
        blob_mask (ndarray): Boolean mask of segmented regions.
        labeled_image (ndarray): Labeled regions, relabeled by ascending longitude.
        merged_props_sorted (list): List of dicts, each containing merged regionprops and lat/lon info.
    """
    from skimage.measure import label, regionprops
    import numpy as np

    # Step 1: Segment image and get regionprops
    blob_mask, labeled_image, props_data, props_intensity = find_blob(
        image_to_segment, intensity_image, threshold_abs, mode
    )

    def rowcol_to_latlon(row, col):
        lat = (90 - lats[0]) - row
        lon = lon_lims[1] - col
        return lat, lon

    # Step 2: Merge each region from props_data and props_intensity by label
    props_by_label = {}
    for region_data in props_data:
        props_by_label[region_data.label] = {
            'label': region_data.label,
            'area': region_data.area,
            'bbox': region_data.bbox,
            'eccentricity': region_data.eccentricity,
            'seg_intensity_max': region_data.intensity_max,
            'seg_intensity_mean': region_data.intensity_mean,
            'seg_intensity_min': region_data.intensity_min#,
        }

        # Coordinates from segmentation image
        if hasattr(region_data, 'centroid'):
            r, c = region_data.centroid
            lat, lon = rowcol_to_latlon(r, c)
            props_by_label[region_data.label]['centroid'] = (r, c)
            props_by_label[region_data.label]['centroid_latlon'] = (lat, lon)
            props_by_label[region_data.label]['centroid_lon'] = lon

        if hasattr(region_data, 'weighted_centroid'):
            r, c = region_data.weighted_centroid
            lat, lon = rowcol_to_latlon(r, c)
            props_by_label[region_data.label]['weighted_centroid'] = (r, c)
            props_by_label[region_data.label]['weighted_centroid_latlon'] = (lat, lon)

        if hasattr(region_data, 'coords'):
            coords_latlon = [rowcol_to_latlon(r, c) for r, c in region_data.coords]
            props_by_label[region_data.label]['coords'] = region_data.coords
            props_by_label[region_data.label]['coords_latlon'] = coords_latlon

    # Step 3: Add intensity image properties
    for region_int in props_intensity:
        label_ = region_int.label
        if label_ in props_by_label:
            props_by_label[label_]['intensity_max'] = region_int.intensity_max
            props_by_label[label_]['intensity_mean'] = region_int.intensity_mean
            props_by_label[label_]['intensity_min'] = region_int.intensity_min

        if hasattr(region_data, 'weighted_centroid'):
            r, c = region_data.weighted_centroid
            lat, lon = rowcol_to_latlon(r, c)
            props_by_label[region_data.label]['weighted_centroid'] = (r, c)
            props_by_label[region_data.label]['weighted_centroid_latlon'] = (lat, lon)

    # Step 4: Sort by longitude
    merged_props = list(props_by_label.values())
    merged_props_sorted = sorted(merged_props, key=lambda x: x['centroid_lon'])

    # Step 5: Relabel regions in labeled image
    new_labeled_image = np.zeros_like(labeled_image)
    label_mapping = {}
    for new_label, region in enumerate(merged_props_sorted, start=1):
        old_label = region['label']
        new_labeled_image[labeled_image == old_label] = new_label
        label_mapping[old_label] = new_label
        region['label'] = new_label  # update to new label

    return blob_mask, new_labeled_image, merged_props_sorted


import csv

def export_regions_to_csv(merged_props_sorted, filepath):
    """
    Export merged region properties to a CSV using only Python's built-in csv module.

    Parameters:
        merged_props_sorted (list of dicts): Output from process_blob.
        filepath (str): Path to write the CSV file.
    """

    # Define desired fields (skip std + duplicates, fix centroid keys, add weighted centroids)
    fields = [
        'label',
        'area',
        'eccentricity',
        'seg_intensity_max',
        'seg_intensity_mean',
        'seg_intensity_min',
        'intensity_max',
        'intensity_mean',
        'intensity_min',
        'centroid_lat',
        'centroid_lon',
        'weighted_centroid_lat_seg',
        'weighted_centroid_lon_seg'
    ]

    # Open CSV file and write header + rows
    with open(filepath, mode='w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fields)
        writer.writeheader()

        for region in merged_props_sorted:
            row = {}

            # Direct entries
            for key in fields:
                if key in region:
                    row[key] = region[key]

            # Derive weighted centroids if needed
            # Segmentation image
            if 'weighted_centroid_latlon' in region:
                row['weighted_centroid_lat_seg'] = region['weighted_centroid_latlon'][0]
                row['weighted_centroid_lon_seg'] = region['weighted_centroid_latlon'][1]

            # Intensity image
            if 'weighted_centroid_latlon_intensity' in region:
                row['weighted_centroid_lat_int'] = region['weighted_centroid_latlon_intensity'][0]
                row['weighted_centroid_lon_int'] = region['weighted_centroid_latlon_intensity'][1]

            # Safely handle centroid split
            if 'centroid_latlon' in region:
                row['centroid_lat'] = region['centroid_latlon'][0]
                row['centroid_lon'] = region['centroid_latlon'][1]

            writer.writerow(row)



import matplotlib.pyplot as plt

def plot_regions_on_axis(
    ax,
    labeled_image,
    merged_props,
    plot_labels=True,
    plot_contours=True,
    plot_masks=False,
    contour_color='white',
    mask_alpha=0.3,
    lats=None,
    lon_lims=None,
):
    import numpy as np
    import matplotlib.pyplot as plt
    from skimage.measure import find_contours

    if lats is None or lon_lims is None:
        raise ValueError("Both 'lats' and 'lon_lims' must be provided to convert to lat-lon coordinates.")

    unique_labels = np.unique(labeled_image)
    unique_labels = unique_labels[unique_labels != 0]  # skip background

    def rowcol_to_latlon(row, col):
        lat = (90 - lats[0]) - row
        lon = lon_lims[1] - col
        return lat, lon

    # Plot shaded masks in lat-lon
    if plot_masks:
        for label in unique_labels:
            mask = labeled_image == label
            rgb = plt.matplotlib.colors.to_rgb(contour_color)
            rgba = (*rgb, mask_alpha)

            rows, cols = np.where(mask)
            for r, c in zip(rows, cols):
                lat, lon = rowcol_to_latlon(r, c)
                ax.add_patch(plt.Rectangle(
                    (lon - 0.5, lat - 0.5), 1.0, 1.0,
                    facecolor=rgba,
                    edgecolor='none',
                    linewidth=0,
                    zorder=2
                ))

    # Plot contours in lat-lon
    if plot_contours:
        for label in unique_labels:
            mask = labeled_image == label
            contours = find_contours(mask.astype(float), 0.5)
            for contour in contours:
                latlon_contour = np.array([rowcol_to_latlon(r, c) for r, c in contour])
                ax.plot(
                    latlon_contour[:, 1],  # longitude
                    latlon_contour[:, 0],  # latitude
                    color=contour_color,
                    linewidth=1.2,
                    alpha=0.8,
                    zorder=3
                )

    # Plot region number labels in lat-lon
    if plot_labels:
        for region in merged_props:
            label = region['label']
            latlon = region.get('centroid_latlon', None)
            if latlon is not None:
                lat, lon = latlon
                ax.text(
                    lon, lat,
                    str(label),
                    color='white',
                    fontsize=8,
                    ha='center',
                    va='center',
                    bbox=dict(boxstyle='round,pad=0.2', fc=contour_color, ec='none', alpha=0.5),
                    zorder=4
                )



def plot_extrema_on_axisa(ax, extrema_dict, data_type, extrema_type, text_color='red',fontsize=8):
    """
    Annotate extrema points on a matplotlib axis with custom text characters.

    Parameters:
    - ax : matplotlib.axes.Axes
        The axis to plot annotations on.
    - extrema_dict : dict
        Output from find_extrema function.
    - data_type : str
        'fNH3', 'PCld', or 'RGB' (or others as added).
    - extrema_type : str
        'max' or 'min'.
    - text_color : str
        Color of text annotations.
    """

    if data_type not in extrema_dict:
        raise ValueError(f"Data type '{data_type}' not found in extrema dictionary.")
    if extrema_type not in ['maxima', 'minima']:
        raise ValueError("extrema_type must be 'maxima' or 'minima'.")

    label_map = {
        ('NH3', 'maxima'): 'N',
        ('NH3', 'minima'): 'D',
        ('PCloud', 'maxima'): 'L',
        ('PCloud', 'minima'): 'H',
        ('RGB',  'maxima'): 'P',
        ('RGB',  'minima'): '5'
    }

    marker_label = label_map.get((data_type, extrema_type), '?')
    extrema_data = extrema_dict[data_type][extrema_type]
    coords = extrema_data['coords']

    if coords is None or len(coords) == 0:
        print(f"No {extrema_type} coordinates to plot for {data_type}.")
        return

    #lats, lons = zip(*coords)

    for lat, lon in coords:
        ax.text(lon, lat, marker_label, color=text_color, 
                horizontalalignment='center', verticalalignment='center', 
                fontsize=fontsize)

    #ax.set_xlabel("Longitude")
    #ax.set_ylabel("Latitude")

def extrema_overplot_all(results,axes = {'axNH3': False, 'axCH4': False, 'axRGB': False}):
    # Define the plotting parameters for each (data_type, extrema_type)
    plot_specs = {
        ('NH3', 'minima'): {'color': 'k', 's': 8, 'lw': 1.0},
        ('NH3', 'maxima'): {'color': 'w', 's': 8, 'lw': 1.0},
        ('PCloud', 'minima'): {'color': 'b', 's': 8, 'lw': 0.5},
        ('PCloud', 'maxima'): {'color': 'y', 's': 8, 'lw': 0.5},
        ('RGB', 'minima'): {'color': 'C1', 's': 8, 'lw': 0.5},
        ('RGB', 'maxima'): {'color': 'C0', 's': 8, 'lw': 0.5},
    }
    
    # Adjust sizes/linewidths by axis if needed
    axis_adjustments = {
        'axNH3': {
            ('NH3', 'minima'): {'s': 10, 'lw': 0.5},
            ('NH3', 'maxima'): {'s': 10, 'lw': 0.5},
        },
        'axCH4': {
            ('PCloud', 'minima'): {'s': 10, 'lw': 1.0},
            ('PCloud', 'maxima'): {'s': 10, 'lw': 1.0},
        },
        'axRGB': {
            ('NH3', 'minima'): {'s': 0, 'lw': 0.5},
            ('RGB', 'minima'): {'s': 0, 'lw': 1.0},
            ('RGB', 'maxima'): {'s': 0, 'lw': 1.0},
        }
    }
    
    # Loop through axes and plot
    
    
    for ax_name, ax in axes.items():
        for (data_type, extrema_type), base_spec in plot_specs.items():
            spec = base_spec.copy()
            # Override with axis-specific adjustments if available
            overrides = axis_adjustments.get(ax_name, {}).get((data_type, extrema_type))
            if overrides:
                spec.update(overrides)
    
            plot_extrema_on_axisa(
                ax, results,
                data_type=data_type,
                extrema_type=extrema_type,
                text_color=spec['color'],
                fontsize=spec['s']#,
                #linewidth=spec['lw']
            )
