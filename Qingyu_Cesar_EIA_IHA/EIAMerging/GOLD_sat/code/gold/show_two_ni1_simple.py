#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 16:56:09 2024

@author: qingyuzhu
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np

# NH
fname='/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/050/GOLD_L1C_CHA_NI1_2023_051_00_22_v05_r01_c01.nc'

# Open in read-only mode
nc = Dataset(fname, 'r')

# Get a list of all variables
print(nc.variables.keys())



wavelengths = nc.variables['WAVELENGTH'][:]

nx,ny,nw=wavelengths.shape

wavelength=np.zeros(nw)
for i in range(nw):
    
    wavelength[i]=np.nanmean(wavelengths[:,:,i])

idx_1356 = np.argmin(np.abs(wavelength - 135.6))


grid_ew = nc.variables['REFERENCE_POINT_LON'][:] * 1
grid_ns = nc.variables['REFERENCE_POINT_LAT'][:] * 1

radiance = nc.variables['RADIANCE'][:]
em_1356 = radiance[:,:, idx_1356]


sza=nc.variables['SOLAR_ZENITH_ANGLE'][:]

em_1356[sza<=95]=np.nan

# Close the file when done
nc.close()

valid_rows = ~np.isnan(grid_ew).all(axis=1)
valid_cols = ~np.isnan(grid_ew).all(axis=0)

# Subset your arrays if that makes sense in your data model:
grid_ew_sub = grid_ew[valid_rows][:, valid_cols]
grid_ns_sub = grid_ns[valid_rows][:, valid_cols]
em_1356_sub = em_1356[valid_rows][:, valid_cols]

# Now these sub-arrays (hopefully) contain no NaNs in x or y
plt.figure()
plt.pcolormesh(grid_ew_sub, grid_ns_sub, em_1356_sub, shading='auto', vmin=0, vmax=600)


#%% SH
fname='/Users/qingyuzhu/Downloads/tmp 2/archive_L1C/2023/051/GOLD_L1C_CHB_NI1_2023_052_00_25_v05_r01_c01.nc'

# Open in read-only mode
nc = Dataset(fname, 'r')

# Get a list of all variables
print(nc.variables.keys())



wavelengths = nc.variables['WAVELENGTH'][:]

nx,ny,nw=wavelengths.shape

wavelength=np.zeros(nw)
for i in range(nw):
    
    wavelength[i]=np.nanmean(wavelengths[:,:,i])

idx_1356 = np.argmin(np.abs(wavelength - 135.6))


grid_ew = nc.variables['REFERENCE_POINT_LON'][:] * 1
grid_ns = nc.variables['REFERENCE_POINT_LAT'][:] * 1

radiance = nc.variables['RADIANCE'][:]
em_1356 = radiance[:,:, idx_1356]


sza=nc.variables['SOLAR_ZENITH_ANGLE'][:]

em_1356[sza<=95]=np.nan

# Close the file when done
nc.close()

valid_rows = ~np.isnan(grid_ew).all(axis=1)
valid_cols = ~np.isnan(grid_ew).all(axis=0)

# Subset your arrays if that makes sense in your data model:
grid_ew_sub = grid_ew[valid_rows][:, valid_cols]
grid_ns_sub = grid_ns[valid_rows][:, valid_cols]
em_1356_sub = em_1356[valid_rows][:, valid_cols]

# Now these sub-arrays (hopefully) contain no NaNs in x or y

plt.pcolormesh(grid_ew_sub, grid_ns_sub, em_1356_sub, shading='auto', vmin=0, vmax=600)


plt.xlim(-80,0)
plt.ylim(-40,40)
plt.colorbar()

plt.show()

