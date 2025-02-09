#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 11:11:28 2025

@author: qingyuzhu
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse

def visualize_wind_data(output_file):
    """Visualize the ut_mid and wind data from an npz file."""
    try:
        # Load the data from the .npz file
        data = np.load(output_file)

        # Extract ut_mid and wind values
        ut_mid1 = data['ut_mid1']
        wind1 = data['wind1']
        
        ut_mid2 = data['ut_mid2']
        wind2 = data['wind2']

        # Plot the data
        plt.figure(figsize=(10, 6))
        plt.plot(ut_mid1, wind1, marker='o', markerfacecolor="None", linestyle='-',  label='Meridional')
        plt.plot(ut_mid2, wind2, marker='o', markerfacecolor="None", linestyle='-',  label='Zonal')
        im=plt.scatter(ut_mid1, wind1, s=16, c= data['qcode1'],vmin=0,vmax=2,cmap=plt.cm.bwr)
        im=plt.scatter(ut_mid2, wind2, s=16, c= data['qcode2'],vmin=0,vmax=2,cmap=plt.cm.bwr)
        plt.colorbar(im)
        plt.xlabel('UT (hours)')
        plt.ylabel('Wind (m/s)')
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()

        # Show the plot
        plt.tight_layout()
        plt.show()

    except FileNotFoundError:
        print(f"Error: File {output_file} not found.")
    except KeyError as e:
        print(f"Error: Missing key in the npz file: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize wind data from an npz file.")
    parser.add_argument(
        "output_file", type=str, help="The full path of the npz file to visualize."
    )
    args = parser.parse_args()

    # Visualize the wind data
    visualize_wind_data(args.output_file)
