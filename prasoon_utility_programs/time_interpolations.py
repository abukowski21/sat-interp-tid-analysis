"""
This will do the time interpolations of electron density
created by Prasoon on 16 September 2023
"""


import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import os, glob
import pandas as pd
from datetime import datetime, timedelta
from scipy.interpolate import CubicSpline
import importlib
from tqdm import tqdm
import pyIGRF as pyigrf
import apexpy
from apexpy import Apex



def edens_sattime(sat_data, start = None, end = None):
    """
    The purpose of the function is to do time interpolations such that 
    from the multiple Edens plots obtained for different Sami times
    we can get a single Edens plots that evolves with sat_step
    """
    
    day_s = sat_data['sat_time'][start].dt.day.values
    hour_s = sat_data['sat_time'][start].dt.hour.values
    day_e = sat_data['sat_time'][end].dt.day.values
    hour_e = sat_data['sat_time'][end].dt.hour.values
    sat_data1 = sat_data.isel(sat_step=slice(start, end+1))
    sat_interval = []

    # Storing all sat_time in seconds for the index interval (start, end+1) of sat_time
    for i in tqdm(sat_data1['sat_time']):
        
        sat_interval.append(i.dt.second.values + 60*i.dt.minute.values + 3600*i.dt.hour.values)

    sami_s = 0
    sami_e = 0
    
    # Finding the sami_time interval corresponding to sat_interval 
    for i in tqdm(range(len(sat_data1['sami_time'])-1)):
        
        date_sami = sat_data1['sami_time'][i]

        # Finding index of sami_time corresponding to the start of sat_interval
        if date_sami.dt.day.values == day_s:
            
            if date_sami.dt.hour.values == hour_s:
                time_sami = date_sami.dt.second.values + 60*date_sami.dt.minute.values + date_sami.dt.hour.values*3600
                
                if (time_sami <= sat_interval[0] and (time_sami+5*60) > sat_interval[0]) or (time_sami >= sat_interval[0] and (time_sami+5*60) < sat_interval[0]):
                    sami_s = i

        # Finding index of sami_time corresponding to the end of sat_interval
        if date_sami.dt.day.values == day_e:
            
            if date_sami.dt.hour.values == hour_e:
                time_sami = date_sami.dt.second.values + date_sami.dt.minute.values*60 + date_sami.dt.hour.values*3600
                
                if (time_sami <= sat_interval[-1] and (time_sami+5*60) > sat_interval[-1]) or (time_sami >= sat_interval[-1] and (time_sami+5*60) < sat_interval[-1]):
                    sami_e = i
                    break

    # Defining a dataset that contains data limited to the satellite pass-only
    sat_data2 = sat_data1.isel(sami_time=slice(sami_s, sami_e+1))
    print(sat_data2)
    
    sami_time_dt = sat_data2['sami_time']
    sami_time = []
    
    for i in tqdm(sami_time_dt):
        
        sami_time.append(i.dt.second.values + i.dt.minute.values*60 + 3600*i.dt.hour.values)
    
    interp_edens = []

    # Applying cubic spline interpolations function
    # Using minutes to be the x-axis and Edens to be the y-axis
    for i in tqdm(range(len(sat_data2['sat_time']))):
        
        edens = []
        
        for j in range(len(sat_data2['sami_time'])):            
            
            edens.append(sat_data2['edens'][j,i].values)
                    
        cs = CubicSpline(sami_time, edens)
        time_i = sat_data2['sat_time'][i]
        
        # Storing values of time in seconds to find values of edens
        time = (time_i.dt.second.values) + 60*(time_i.dt.minute.values) + 3600*(time_i.dt.hour.values)
        interp_edens.append(cs(time))

    # Creating new data variable storing edens 
    sat_data2['interp_edens'] = ('sat_step', interp_edens)
    return sat_data2

    
def magnetic_coords(sat_data):
    mlat = []
    mlon = []
    
    # Introducing the MLAT, MLON, MLT columns in satlocddf
    for glat, glon, alt, gtime in tqdm(zip(sat_data.glat.values, sat_data.glon.values, sat_data.alt.values, sat_data.sat_time.values), total=len(sat_data.glon)):
        
        date_str = (np.datetime_as_string(gtime))[:19]
        format_string = "%Y-%m-%dT%H:%M:%S"
        dtime = datetime.strptime(date_str, format_string)
        decimal_year = dtime.year + ((dtime - datetime(dtime.year, 1, 1)).days)/365

        apex = Apex(decimal_year)
        mlatitude, mlongitude = apex.convert(glat, glon, 'geo', 'apex', height = alt)        
        #mlatitude, mlongitude = magfield
        mlat.append(mlatitude)
        mlon.append(mlongitude)

    # Adding mlat and mlon data variables in dataset
    sat_data['mlat'] = ('sat_step', mlat)
    sat_data['mlon'] = ('sat_step', mlon)
    
    return sat_data


def plot_edens(sat1, sat2, sat3, sat4, sat5, sat6, variable):

    # Plotting the edens on single plot. Variable is a string defining the data variable being plotted against mlat
    c = ['black', 'purple', 'blue', 'yellow', 'orange', 'red']
    label = ['quiet 00:20 UTC', 'quiet 02:02 UTC', 'quiet 03:44 UTC', 'storm 15:39 UTC', 'storm 17:19 UTC', 'storm 19:01 UTC']
    
    plt.figure(figsize=(12,8))
    plt.plot(sat1['mlat'], sat1[variable], color = c[0], label = label[0])
    plt.plot(sat2['mlat'], sat2[variable], color = c[1], label = label[1])
    plt.plot(sat3['mlat'], sat3[variable], color = c[2], label = label[2])
    plt.plot(sat4['mlat'], sat4[variable], color = c[3], label = label[3])
    plt.plot(sat5['mlat'], sat5[variable], color = c[4], label = label[4])
    plt.plot(sat6['mlat'], sat6[variable], color = c[5], label = label[5])
    plt.legend(title='Satellite passes', loc = 'upper left')
    plt.xlabel('Magnetic Latitude MLAT')
    plt.ylabel(variable)
    plt.title('Comparison Plots for Quiet and Storm Time Satellite Passes')
    plt.show()
    
    return
