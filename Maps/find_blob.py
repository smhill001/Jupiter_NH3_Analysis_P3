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


def convert_one_time(row, col, mean_time_array):
    #print("%%%%%%%%%%row_inds, col_inds=",row, col)
    #print("mean_time_array[row_inds[0], col_inds[0]]=",Time(mean_time_array[row, col], format='jd'))
    try:
        #print("mean_time_array[row, col]=",mean_time_array[row, col],row,col)
        #temptime = Time(mean_time_array[row, col], format='jd')
        temptime = Time(mean_time_array[row, col], format='jd')
        fits_time=temptime.fits
        tempJD = deepcopy(temptime)
        tempJD.format = 'jd'
        jd_time=tempJD.value
    except:
        fits_time='N/A'
        jd_time=0.0
    return fits_time, jd_time


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


def process_blob(image_to_segment, intensity_image, lats, lon_lims, timearray=None, threshold_abs=None, mode='max'):
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
            
        # Times from timearray
        if timearray is not None:
            # Time info
            #row_inds = coords_sorted[:, 0]
            #col_inds = coords_sorted[:, 1]
            ###!!!!!!! Need to get these coords sorted out and maybe trim the
            ###!!!!!!! mean-time array to the proper patch or simply to 
            ###!!!!!!! a 1D longitude array. The "15" below should not be a
            ###!!!!!!! fixed value!
            #templat,templon=15-coords_sorted[:, 0],coords_sorted[:, 1]
            #row_inds, col_inds=90-templat,coords_sorted[:, 1]


            r, c = region_data.weighted_centroid
            row, col = int(round(r)), int(round(c))
            fits_time, jd_time = convert_one_time(row, col, timearray)

            #print("###########",row, col)
            #print("###########",fits_time, jd_time)
            try:
                #jd_value = timearray[wy, wx]
                #time_obj = Time(jd_value, format='jd')
                #fits_time = time_obj.to_value('isot', subfmt='date_hms')
                props_by_label[region_data.label]['times'] = {
                    'jd_time': jd_time,
                    'fits_time': fits_time
                }
            except Exception as e:
                props_by_label[region_data.label]['times'] = {
                    'jd_time': None,
                    'fits_time': None,
                    'error': str(e)
                }

        #merged_props.append(record)


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
        'jd_time',
        'fits_time',
        'area',
        'eccentricity',
        'centroid_lat',
        'centroid_lon',
        'weighted_centroid_lat_seg',
        'weighted_centroid_lon_seg',
        'seg_intensity_mean',
        'seg_intensity_min',
        'seg_intensity_max',
        'intensity_mean',
        'intensity_min',
        'intensity_max'
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

            if 'jd_time' in region['times']:
                row['jd_time'] = region['times']['jd_time']
            if 'fits_time' in region['times']:
                row['fits_time'] = region['times']['fits_time']


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
