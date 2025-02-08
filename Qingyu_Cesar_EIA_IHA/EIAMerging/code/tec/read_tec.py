#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 08:33:30 2024

@author: qingyuzhu

Read TEC
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# File path to your data
file_path = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/tec/gps230412g.002.txt'

# Define the column names
columns = [
    "YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "RECNO", "KINDAT", "KINST", 
    "UT1_UNIX", "UT2_UNIX", "GDLAT", "GLON", "TEC", "DTEC"
]

# Read the file, skipping the first row since it contains text
data = pd.read_csv(file_path, sep=r'\s+', names=columns, skiprows=1)

# Combine year, month, day, hour, and minute into a single datetime column
# Rename columns to expected format
data = data.rename(columns={"YEAR": "year", "MONTH": "month", "DAY": "day", "HOUR": "hour", "MIN": "minute"})

# Combine the columns into a single datetime column
data['DATETIME'] = pd.to_datetime(data[["year", "month", "day", "hour", "minute"]])

# Get unique datetimes
unique_datetimes = data['DATETIME'].unique()

#%%

# Loop through each unique datetime
for datetime in unique_datetimes:
    

    # Filter the data for the current datetime
    subset = data[data['DATETIME'] == datetime]
    
    # Create a pivot table for TEC as a function of GDLAT and GLON
    pivot = subset.pivot(index='GDLAT', columns='GLON', values='TEC')
    
    # Skip if the pivot table is empty
    if pivot.empty:
        continue
    
    # Plotting
    plt.figure(figsize=(10, 6))
    plt.pcolormesh(
        pivot.columns, pivot.index, pivot.values, 
        shading='auto', cmap='viridis', vmin=0, vmax=60)
    
    plt.colorbar(label='TEC (Total Electron Content)')
    plt.title(f"TEC as a function of GLON and GDLAT\nDatetime: {datetime}")
    plt.xlabel("GLON (Longitude)")
    plt.ylabel("GDLAT (Latitude)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    
    

