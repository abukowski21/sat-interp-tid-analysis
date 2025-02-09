#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 28 07:23:44 2024

@author: qingyuzhu
"""

from datetime import datetime
import re
import pandas as pd
import os

# Extract timestamps from filenames
def extract_time(filename):
    match = re.search(r"_(\d{4})_(\d{3})_(\d{2})_(\d{2})_", filename)
    if match:
        year, doy, hour, minute = map(int, match.groups())
        timestamp = datetime.strptime(f"{year} {doy} {hour} {minute}", "%Y %j %H %M")
        return timestamp
    return None

# Pair CHA and CHB filenames by closest timestamps
def find_closest_pairs(cha_files, chb_files):
    cha_times = [(filename, extract_time(filename)) for filename in cha_files]
    chb_times = [(filename, extract_time(filename)) for filename in chb_files]
    closest_pairs = []

    for cha_file, cha_time in cha_times:
        closest_chb = min(chb_times, key=lambda x: abs(x[1] - cha_time))
        closest_pairs.append((cha_file, closest_chb[0], abs(closest_chb[1] - cha_time)))

    return closest_pairs

#%%
rewrite_file = -1 

# Base directory containing year and day-of-year folders
base_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/'

# Directory to save output CSV files
output_dir = '/Volumes/Io2/proposal/2025/PostSunsetEIA/data/L1C/pairs/'
os.makedirs(output_dir, exist_ok=True)

# Iterate over years and days of year
for year in range(2022, 2023):  # Adjust range as needed
    for doy in range(1, 366):  # 1 to 365 (or 366 for leap years)
        data_dir = os.path.join(base_dir, f"{year}/{doy:03d}/")
        if not os.path.exists(data_dir):
            continue
        
        output_filename = os.path.join(output_dir, f"closest_pairs_{year}_{doy:03d}.csv")
        if (os.path.exists(output_filename)) & (rewrite_file<0):
            print(f"Processed file {output_filename} already exists. Skipping...")
            continue

        # List CHA and CHB files
        cha_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if "CHA" in f and f.endswith(".nc") and not f.startswith("._")]
        chb_files = [os.path.join(data_dir, f) for f in os.listdir(data_dir) if "CHB" in f and f.endswith(".nc") and not f.startswith("._")]

        if not cha_files or not chb_files:
            continue

        # Find closest pairs
        closest_pairs = find_closest_pairs(cha_files, chb_files)

        # Create a DataFrame for the day's results
        day_results_df = pd.DataFrame(closest_pairs, columns=["CHA File", "CHB File", "Time Difference"])
        day_results_df.sort_values("Time Difference", inplace=True)

        # Save results to a CSV file for the specific day
        
        day_results_df.to_csv(output_filename, index=False)

        print(f"Closest pairs for {year} DOY {doy} saved to '{output_filename}'")

