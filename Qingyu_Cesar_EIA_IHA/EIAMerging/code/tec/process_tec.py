#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 09:23:36 2024

@author: qingyuzhu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob

rewrite_file = -1

# Directory containing the data files
data_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/tec/temp/'
csv_output_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/tec/csv/'
plot_output_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/tec/plots/'

# Ensure the output directories exist
os.makedirs(csv_output_dir, exist_ok=True)
os.makedirs(plot_output_dir, exist_ok=True)

# Define the column names
columns = [
    "YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "RECNO", "KINDAT", "KINST", 
    "UT1_UNIX", "UT2_UNIX", "GDLAT", "GLON", "TEC", "DTEC"
]

# List all files in the folder
file_list = sorted(glob(os.path.join(data_dir, '*.txt')))

# Loop through each file
for file_path in file_list:
    
    print ()
    print(f"Processing file: {file_path}")
    
    
    date_key='20'+os.path.basename(file_path)[3:9]
    
    

    csv_output_dir1 = os.path.join(csv_output_dir, date_key)
    plot_output_dir1 = os.path.join(plot_output_dir, date_key)
    
    
    
    
    if (os.path.exists(csv_output_dir1))& (rewrite_file<0):
        print(f"Processed file {file_path} already exists. Skipping...")
        continue
    
    os.makedirs(csv_output_dir1, exist_ok=True)
    os.makedirs(plot_output_dir1, exist_ok=True)
    
    # Read the file, skipping the first row since it contains text
    data = pd.read_csv(file_path, sep=r'\s+', names=columns, skiprows=1)

    
    # Combine year, month, day, hour, and minute into a single datetime column
    data = data.rename(columns={"YEAR": "year", "MONTH": "month", "DAY": "day", "HOUR": "hour", "MIN": "minute"})
    data['DATETIME'] = pd.to_datetime(data[["year", "month", "day", "hour", "minute"]])
    
    # Get unique datetimes
    unique_datetimes = data['DATETIME'].unique()
    
    # Loop through each unique datetime
    for datetime in unique_datetimes:
        print(f"Processing datetime: {datetime}")
        
        csv_filename = os.path.join(csv_output_dir1, f"{datetime.strftime('%Y%m%d_%H%M')}.csv")
        plot_filename = os.path.join(plot_output_dir1, f"{datetime.strftime('%Y%m%d_%H%M')}.png")
        
        
        
        # Filter the data for the current datetime
        subset = data[data['DATETIME'] == datetime]
        
        # Create a pivot table for TEC as a function of GDLAT and GLON
        pivot = subset.pivot(index='GDLAT', columns='GLON', values='TEC')
        
        # Skip if the pivot table is empty
        if pivot.empty:
            continue
        
        # Save the data as CSV
        
        subset[['GDLAT', 'GLON', 'TEC']].to_csv(csv_filename, index=False)
        print(f"Saved data to {csv_filename}")
        
        # Plotting
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(
            pivot.columns, pivot.index, pivot.values, 
            shading='auto', cmap='viridis', vmin=0, vmax=60
        )
        plt.colorbar(label='TEC (Total Electron Content)')
        plt.title(f"TEC as a function of GLON and GDLAT\nDatetime: {datetime}")
        plt.xlabel("GLON (Longitude)")
        plt.ylabel("GDLAT (Latitude)")
        plt.grid(True)
        plt.tight_layout()
        
        # Save the plot
        
        plt.savefig(plot_filename)
        print(f"Saved plot to {plot_filename}")
        plt.close()
