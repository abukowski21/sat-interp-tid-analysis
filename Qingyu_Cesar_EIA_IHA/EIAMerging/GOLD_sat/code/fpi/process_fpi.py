#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 21:25:33 2024

@author: qingyuzhu
"""

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
import os
import datetime as dt
import pandas as pd

rewrite_file = -1 

# List of stations
stations = ['aif','bfp', 'clf']

# Base path for data
base_path = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/fpi/'
save_base_path = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/fpi/processed/'
save_fig_path = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/fpi/fig/'


for station in stations:
    print ()
    print(f"Processing data for station: {station}")

    fpath = os.path.join(base_path, station)
    save_dir = os.path.join(save_base_path, station)
    save_dir1 = os.path.join(save_fig_path, station)
    os.makedirs(save_dir, exist_ok=True)  # Ensure the save directory exists
    os.makedirs(save_dir1, exist_ok=True)  # Ensure the save directory exists

    flist = sorted(glob(fpath + '/*.txt'))

    for ifile, fname in enumerate(flist):
        print(f"Processing file: {fname}")
        
        name = fname.split('/')[-1][:9]
        save_path = os.path.join(save_dir, f"{name}.csv")
        save_path1 = os.path.join(save_dir1, f"{name}.png")
        
        if (os.path.exists(save_path)) & (rewrite_file<0):
            print(f"Processed file {save_path} already exists. Skipping...")
            continue

        columns = [
            "YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "RECNO", "KINDAT", "KINST",
            "UT1_UNIX", "UT2_UNIX", "ALTB", "ALTE", "GDLAT", "GLON",
            "VN1", "VN2", "DVN1", "DVN2", "GVN1", "GVN2", "DGVN1", "DGVN2", "FPI_DATAQUAL"
        ]
        df = pd.read_csv(fname, sep='\s+', names=columns, skiprows=1)
        
        # Select relevant columns
        df = df[["YEAR", "MONTH", "DAY", "HOUR", "MIN", "VN1", "DVN1", "VN2", "DVN2"]]
        
        # Combine year, month, day, hour, and minute into a single datetime column
        # Rename columns to expected format
        df = df.rename(columns={"YEAR": "year", "MONTH": "month", "DAY": "day", "HOUR": "hour", "MIN": "minute"})
        
        # Combine the columns into a single datetime column
        df['DATETIME'] = pd.to_datetime(df[["year", "month", "day", "hour", "minute"]])
        
        
        # Drop any rows with NaNs in VN1 or VN2
        #df = df.dropna(subset=["VN1", "VN2"])
        
        # Plot the time series
        plt.figure(figsize=(10, 6))
        
        # VN1 with error bars
        plt.errorbar(
            df['DATETIME'], df['VN1'], yerr=df['DVN1'], fmt='o', label='VN1', color='blue'
        )
        
        # VN2 with error bars
        plt.errorbar(
            df['DATETIME'], df['VN2'], yerr=df['DVN2'], fmt='o', label='VN2', color='orange'
        )
        
        # Apply smoothing for VN1 (time-based window)
        df = df.set_index('DATETIME')  # Set the datetime column as the index for rolling()
        df['VN1_SMOOTH'] = df['VN1'].rolling('40min', center=True).mean()  # 20-minute rolling window
        df['VN2_SMOOTH'] = df['VN2'].rolling('40min', center=True).mean()  # 20-minute rolling window
        
        # Reset the index for plotting
        df = df.reset_index()
        


        # Save processed data to CSV
        df.to_csv(save_path, index=False)
        print(f"Processed data saved to {save_path}")
        
        # Plot smoothed VN1
        plt.plot(df['DATETIME'], df['VN1_SMOOTH'], color='blue', label='Zonal', linewidth=4)
        plt.plot(df['DATETIME'], df['VN2_SMOOTH'], color='orange', label='Meridional', linewidth=4)

        # Customize plot
        plt.ylim([-200, 200])
        plt.xlim(df['DATETIME'].min(), df['DATETIME'].max())
        plt.xlabel("Datetime")
        plt.ylabel("Velocity (VN)")
        plt.title(f"({station.upper()}: {name})")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        # Show plot
        plt.savefig(save_path1)
        plt.close()