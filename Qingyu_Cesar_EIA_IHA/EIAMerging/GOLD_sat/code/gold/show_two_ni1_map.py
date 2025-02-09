#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 17:12:41 2024

@author: qingyuzhu
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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

def main():
    # Paths to your NI1 files
    fname_nh = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/2019/280/GOLD_L1C_CHB_NI1_2019_280_20_40_v05_r01_c01.nc'
    fname_sh = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/2019/280/GOLD_L1C_CHB_NI1_2019_280_20_55_v05_r01_c01.nc'

    # Read & subset the NH data
    grid_ew_nh, grid_ns_nh, em_1356_nh = read_and_subset_ni1(fname_nh, sza_threshold=95)

    # Read & subset the SH data
    grid_ew_sh, grid_ns_sh, em_1356_sh = read_and_subset_ni1(fname_sh, sza_threshold=95)

    # Create a Cartopy figure/axis
    fig = plt.figure(figsize=(9, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())

    # Set the map extent (lon_min, lon_max, lat_min, lat_max)
    # Adjust these to focus on the region of interest.
    ax.set_extent([-100, 20, -30, 30], crs=ccrs.PlateCarree())

    # Plot the NH data
    mesh_nh = ax.pcolormesh(grid_ew_nh, grid_ns_nh, em_1356_nh,
                            transform=ccrs.PlateCarree(),
                            shading='auto',
                            vmin=0, vmax=20,
                            cmap='plasma')

    # Plot the SH data
    mesh_sh = ax.pcolormesh(grid_ew_sh, grid_ns_sh, em_1356_sh,
                            transform=ccrs.PlateCarree(),
                            shading='auto',
                            vmin=0, vmax=20,
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
    
    #instrument latitude: 18.3298
    #instrument longitude: 294.6932
    
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
    cbar=fig.colorbar(mesh_sh, ax=ax, orientation='vertical', pad=0.02)
    cbar.set_label('135.6 nm Emission', fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    #plt.title("NASA GOLD NI1 135.6 nm Emission (NH & SH)")
    ax.set_title("NASA GOLD NI1 135.6 nm Emission (NH & SH)",fontsize=16)
    plt.show()

if __name__ == "__main__":
    main()
