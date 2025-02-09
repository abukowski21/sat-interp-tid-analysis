#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 07:46:48 2024

@author: qingyuzhu
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from datetime import datetime
import datetime as dt
import re
import os
import pandas as pd
#%%
def read_and_subset_ni1(fname, sza_threshold=95.0):
    """
    Reads a GOLD L1C NI1 file, extracts 135.6 nm emission (em_1356) and 
    the REFERENCE_POINT_LON/REFERENCE_POINT_LAT grids, masks daylit areas
    (SZA <= sza_threshold), and returns sub-grids with non-NaN columns/rows removed.

    :param fname: Full path to the GOLD L1C NI1 NetCDF file.
    :param sza_threshold: Solar zenith angle cutoff below which we mask out data (daylit).
    :return: (grid_ew_sub, grid_ns_sub, em_1356_sub) 2D arrays of lon, lat, and emission.
    """
    with Dataset(fname, 'r') as nc:
        # Identify index for 135.6 nm
        wavelengths = nc.variables['WAVELENGTH'][:]  # shape (nx, ny, nw)
        nx, ny, nw = wavelengths.shape
        
        # Create a 1D array of wavelengths by averaging each spectral slice
        wv_1d = np.array([
            np.nanmean(wavelengths[:, :, i]) for i in range(nw)
        ])
        idx_1356 = np.argmin(np.abs(wv_1d - 135.6))

        # Read the coordinate grids (lon & lat)
        grid_ew = nc.variables['REFERENCE_POINT_LON'][:]  # shape (nx, ny)
        grid_ns = nc.variables['REFERENCE_POINT_LAT'][:]  # shape (nx, ny)

        # Extract the 135.6 nm radiance
        radiance = nc.variables['RADIANCE'][:]  # shape (nx, ny, nw)
        em_1356 = radiance[:, :, idx_1356]      # shape (nx, ny)

        # Mask out daylit areas
        sza = nc.variables['SOLAR_ZENITH_ANGLE'][:]  # shape (nx, ny)
        em_1356[sza <= sza_threshold] = np.nan

    # Remove rows/columns that are entirely NaN in the grid_ew
    valid_rows = ~np.isnan(grid_ew).all(axis=1)
    valid_cols = ~np.isnan(grid_ew).all(axis=0)

    # Subset each array
    grid_ew_sub  = grid_ew[valid_rows][:, valid_cols]
    grid_ns_sub  = grid_ns[valid_rows][:, valid_cols]
    em_1356_sub  = em_1356[valid_rows][:, valid_cols]

    return grid_ew_sub, grid_ns_sub, em_1356_sub

#%%
def extract_datetime_from_filename(filename):
    """Extract datetime from GOLD L1C filename."""
    match = re.search(r"_(\d{4})_(\d{3})_(\d{2})_(\d{2})_", filename)
    if match:
        year, doy, hour, minute = map(int, match.groups())
        date = datetime.strptime(f"{year} {doy} {hour} {minute}", "%Y %j %H %M")
        return date
    return None

#%%
def plot_ni1_maps(cha_data, chb_data, cha_file, chb_file, output_figure_path):
    """
    Plots the NI1 data for CHA and CHB in one plot.

    Parameters:
    - cha_data: Data from CHA file 
    - chb_data: Data from CHB file 
    - cha_file: File name for CHA data.
    - chb_file: File name for CHB data.
    - output_figure_path: Path to save the output figure.
    """
    # Extract average datetime
    cha_datetime = extract_datetime_from_filename(cha_file)
    chb_datetime = extract_datetime_from_filename(chb_file)
    if cha_datetime and chb_datetime:
        avg_datetime = cha_datetime + (chb_datetime - cha_datetime) / 2
        title = avg_datetime.strftime("%Y-%m-%d %H:%M")
    else:
        title = "CHA and CHB Data Combined"
        
    # Read & subset the NH data
    grid_ew_nh, grid_ns_nh, em_1356_nh = read_and_subset_ni1(cha_data, sza_threshold=95)

    # Read & subset the SH data
    grid_ew_sh, grid_ns_sh, em_1356_sh = read_and_subset_ni1(chb_data, sza_threshold=95)
    
    fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    
    ax.set_title(title,fontsize=16)
    ax.set_extent([-100, 20, -30, 30], crs=ccrs.PlateCarree())

    # Plot the NH data
    mesh_nh = ax.pcolormesh(grid_ew_nh, grid_ns_nh, em_1356_nh,
                            transform=ccrs.PlateCarree(),
                            shading='auto',
                            vmin=0, vmax=80,
                            cmap='plasma')

    # Plot the SH data
    mesh_sh = ax.pcolormesh(grid_ew_sh, grid_ns_sh, em_1356_sh,
                            transform=ccrs.PlateCarree(),
                            shading='auto',
                            vmin=0, vmax=80,
                            cmap='plasma')

    # Add coastlines, borders, land, etc.
    ax.add_feature(cfeature.COASTLINE, linewidth=2.0)
    #ax.add_feature(cfeature.BORDERS,   linestyle=':')
    #ax.add_feature(cfeature.LAND,      facecolor='lightgray')
    ax.add_feature(cfeature.OCEAN,     facecolor='white')
    
    # Add gridlines
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,      # <-- enables latitude/longitude labels
        linewidth=1,
        color='gray',
        alpha=0.5,
        linestyle='--'
    )
    
    # Turn off top/right labels (if desired)
    gl.top_labels    = False
    gl.right_labels  = False
    
    # Format lat/lon ticks nicely
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    # Optionally, customize label style
    gl.xlabel_style = {'size': 16, 'color': 'black'}
    gl.ylabel_style = {'size': 16, 'color': 'black'}

    # Coordinates to highlight
    lat_spot = -22.703764
    lon_spot = 314.99071  # approx. 314.99071° E
    
    # Add a marker (red star) at that spot
    ax.plot(
        lon_spot, lat_spot,
        marker='*',      # or 'o' for a circle, '^' for triangle, etc.
        markersize=15,   # Adjust size to your preference
        color='magenta',
        transform=ccrs.PlateCarree(),
        label='Highlight Spot'
    )
    
    ax.text(
        lon_spot - 15,  # shift a bit in lon so text is not on top of the marker
        lat_spot - 5,
        "Cachoeira Paulista",
        transform=ccrs.PlateCarree(),
        fontsize=16,
        color='magenta',
        bbox=dict( fc="None", ec="None", alpha=0.8)
    )
    
    lat_spot = 18.3298
    lon_spot = 294.6932  # approx. 314.99071° E
    
    # Add a marker (red star) at that spot
    ax.plot(
        lon_spot, lat_spot,
        marker='*',      # or 'o' for a circle, '^' for triangle, etc.
        markersize=15,   # Adjust size to your preference
        color='magenta',
        transform=ccrs.PlateCarree(),
        label='Highlight Spot'
    )
    
    ax.text(
        lon_spot - 5,  # shift a bit in lon so text is not on top of the marker
        lat_spot + 3,
        "Culebra",
        transform=ccrs.PlateCarree(),
        fontsize=16,
        color='magenta',
        bbox=dict( fc="None", ec="None", alpha=0.8)
    )
    
    # Add a colorbar (attaching to one of the mesh objects is fine)
    cbar = fig.colorbar(mesh_sh, ax=ax, orientation='vertical', pad=0.02, aspect=50)
    cbar.set_label('135.6 nm Emission', fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    plt.tight_layout()
    plt.savefig(output_figure_path+avg_datetime.strftime("%Y%m%d_%H%M"))
    plt.close(fig)


#%%
def adjust_doy_for_filename(datetime_obj, threshold_hour=6):
    """
    Adjusts the DOY in the filename based on a time threshold.
    If the hour is less than the threshold, subtract one day.

    :param datetime_obj: Datetime object extracted from the filename.
    :param threshold_hour: Hour threshold to decide the DOY adjustment.
    :return: Adjusted year and DOY as integers.
    """
    if datetime_obj.hour < threshold_hour:
        datetime_obj -= dt.timedelta(days=1)
    return datetime_obj.year, datetime_obj.timetuple().tm_yday
#%%
rewrite_file = -1 

# Directories
input_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/pairs/'  # Directory containing the CSV files
output_fig_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/figures/'  # Directory to save figures
os.makedirs(output_fig_dir, exist_ok=True)

# List all CSV files in the directory
csv_files = [f for f in os.listdir(input_dir) if f.endswith(".csv") and not f.startswith("._")]

for csv_file in csv_files:
    csv_path = os.path.join(input_dir, csv_file)
    df = pd.read_csv(csv_path)

    npair = len(df)
    
    for index, row in df.iterrows():
        cha_file = row["CHA File"]
        chb_file = row["CHB File"]
        
        
        cha_datetime = extract_datetime_from_filename(cha_file)
        chb_datetime = extract_datetime_from_filename(chb_file)

        avg_datetime = cha_datetime + (chb_datetime - cha_datetime) / 2
        
        if avg_datetime.hour<6:
            avg_datetime1=avg_datetime - dt.timedelta(days=1) 
        else:
            avg_datetime1=avg_datetime
            
            
            
        save_year=avg_datetime1.year 
        save_doy =avg_datetime1.timetuple().tm_yday
        

        
        output_figure_path = os.path.join(output_fig_dir, f"{save_year}/{save_doy:03d}/")
        os.makedirs(output_figure_path, exist_ok=True)
        
        fig_files = [f for f in os.listdir(output_figure_path) if f.endswith(".png") and not f.startswith("._")]
        

        if (len(fig_files)==npair)&(rewrite_file<0):
            continue
        
        try:
            plot_ni1_maps(cha_file, chb_file, cha_file, chb_file, output_figure_path)
            print ()
            print(f"Figure created for {cha_file} and {chb_file}, saved to {output_figure_path}")
        except Exception as e:
            print ()
            print(f"Failed to process {cha_file} and {chb_file}: {e}")
