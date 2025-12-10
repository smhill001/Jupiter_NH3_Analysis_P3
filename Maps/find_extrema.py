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
    #print("%%%%%%%%%%row_inds, col_inds=",row_inds, col_inds)
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


def find_extrema(data, min_distance=4, threshold_abs=None, mode='max'):
    if threshold_abs is None:
        threshold_abs = np.nanmean(data)

    if mode == 'min':
        extrema = peak_local_max(-data, min_distance=min_distance, threshold_abs=-threshold_abs)
    else:
        extrema = peak_local_max(data, min_distance=min_distance, threshold_abs=threshold_abs)

    return extrema


def process_extrema(data_arrays: dict, mean_time_array, lats, lon_lims, min_distance=4):
    results = {}

    for label, data in data_arrays.items():
        data = np.array(data)
        data[data == 0] = np.nan
        threshold = np.nanmean(data)
        print(label)
        extrema_result = {}
        for mode in ['max', 'min']:
            extrema = find_extrema(data, min_distance=min_distance, threshold_abs=threshold, mode=mode)

            values = {key: arr[extrema[:, 0], extrema[:, 1]] for key, arr in data_arrays.items()}
            coords_transformed = np.zeros_like(extrema)
            coords_transformed[:, 0] = (90 - lats[0]) - extrema[:, 0]
            coords_transformed[:, 1] = lon_lims[1] - extrema[:, 1]

            # Sort by longitude
            sort_idx = coords_transformed[:, 1].argsort()
            transformed_sorted = coords_transformed[sort_idx]
            coords_sorted = extrema[sort_idx]
            values_sorted = {k: v[sort_idx] for k, v in values.items()}

            # Time info
            row_inds = coords_sorted[:, 0]
            col_inds = coords_sorted[:, 1]
            ###!!!!!!! Need to get these coords sorted out and maybe trim the
            ###!!!!!!! mean-time array to the proper patch or simply to 
            ###!!!!!!! a 1D longitude array. The "15" below should not be a
            ###!!!!!!! fixed value!
            #templat,templon=15-coords_sorted[:, 0],coords_sorted[:, 1]
            #row_inds, col_inds=90-templat,coords_sorted[:, 1]

            fits_times, jd_times = convert_time(row_inds, col_inds, mean_time_array)

            extrema_result[f"{mode}ima"] = {
                "coords": transformed_sorted,
                "pixels": coords_sorted,
                "values": values_sorted,
                "times": {
                    "fits": fits_times,
                    "jd": jd_times
                }
            }

        results[label] = extrema_result

    return results



def export_extrema_to_csv(results, output_filename):
    with open(output_filename, mode='w', newline='') as file:
        writer = None

        for label, extrema_types in results.items():
            for extrema_type, data in extrema_types.items():
                # Compose the header once
                if writer is None:
                    header = [
                        "Label", "Extrema_Type", "Index", "Time_FITS", "Time_JD",
                        "Pixel_Row", "Pixel_Col", "Latitude", "Longitude"
                    ] + [f"Value_{k}" for k in data['values'].keys()]
                    writer = csv.writer(file)
                    writer.writerow(header)

                for i in range(len(data['coords'])):
                    row = [
                        label,
                        extrema_type,
                        i,
                        data['times']['fits'][i],
                        data['times']['jd'][i],
                        data['pixels'][i][0],
                        data['pixels'][i][1],
                        data['coords'][i][0],
                        data['coords'][i][1]
                    ] + [data['values'][k][i] for k in data['values'].keys()]
                    writer.writerow(row)


import matplotlib.pyplot as plt

def plot_extrema_on_axis(ax, extrema_dict, data_type, extrema_type, symbol_color='red',s=15,linewidth=0.5):
    """
    Plot extrema points on a given matplotlib axis.

    Parameters:
    - ax : matplotlib.axes.Axes
        The axis object to plot on.
    - extrema_dict : dict
        Dictionary of extrema results from find_extrema.
    - data_type : str
        'fNH3', 'PCld', or any other data type used in the dictionary.
    - extrema_type : str
        Either 'max' or 'min' to select which extrema to plot.
    - symbol_color : str
        Color of the plotted symbols (matplotlib color string).
    """

    if data_type not in extrema_dict:
        raise ValueError(f"Data type '{data_type}' not found in extrema dictionary.")
    if extrema_type not in ['maxima', 'minima']:
        raise ValueError("extrema_type must be 'maxima' or 'minima'.")

    extrema_data = extrema_dict[data_type][extrema_type]
    coords = extrema_data['coords']

    if coords is None or len(coords) == 0:
        print(f"No {extrema_type} coordinates to plot for {data_type}.")
        return

    # Plotting with '+' for max and '–' for min
    marker = '+' if extrema_type == 'maxima' else '_'

    lats, lons = zip(*coords)
    ax.scatter(lons, lats, marker=marker, color=symbol_color, s=s,linewidth=linewidth,label=f'{data_type} {extrema_type}')

    #ax.set_xlabel("Longitude")
    #ax.set_ylabel("Latitude")
    #ax.legend()

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
    # !!!!!!!NEED FIX FOR ADDITIONAL AXES, E.G., 5um and 889CH4
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
