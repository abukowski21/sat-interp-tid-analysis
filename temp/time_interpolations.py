"""This will do the time interpolations of electron density

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



def edens_sattime(sat_data, sat_start = None, sat_end = None):
    """
    The purpose of function is to do time interpolations such that 
    from the multiple edens plots obtained for different sami times
    we can obtain a single edens plots that evolves with sat_time
    """
    start_day = sat_data['sat_time'][sat_start].dt.day.values
    start_hour = sat_data['sat_time'][sat_start].dt.hour.values
    end_day = sat_data['sat_time'][sat_end].dt.day.values
    end_hour = sat_data['sat_time'][sat_end].dt.hour.values
    sat_data_trim = sat_data.isel(sat_time=slice(sat_start, sat_end+1))
    passing_time = []

    for i in tqdm(sat_data_trim['sat_time']):
        passing_time.append(i.dt.minute.values + 60*i.dt.hour.values)

    sami_start = 0
    sami_end = 0
    
    for i in tqdm(range(len(sat_data_trim['sami_time'])-1)):
        # Finding the sami_time range where sat_time shows pass of satellite 
        date_sami = sat_data_trim['sami_time'][i]
        
        if date_sami.dt.day.values == start_day:
            
            if date_sami.dt.hour.values == start_hour:
                minute_sami = date_sami.dt.minute.values + date_sami.dt.hour.values*60
                if (minute_sami <= passing_time[0] and (minute_sami+5) > passing_time[0]) or (minute_sami >= passing_time[0] and (minute_sami+5) < passing_time[0]):
                    sami_start = i
                    
        if date_sami.dt.day.values == end_day:
            
            if date_sami.dt.hour.values == end_hour:
                minute_sami = date_sami.dt.minute.values + date_sami.dt.hour.values*60
                if (minute_sami <= passing_time[-1] and (minute_sami+5) > passing_time[-1]) or (minute_sami >= passing_time[-1] and (minute_sami+5) < passing_time[-1]):
                    sami_end = i
                    break

    # Defining dataset which contains data limited to the satellite pass only
    sat_final = sat_data_trim.isel(sami_time=slice(sami_start, sami_end+1))
    
    # Applying cubic spline interpolations function
    # Using minutes to be x-axis and edens to be y-axis
    sami_pass = sat_final['sami_time']
    sami_pass_min = []
    for i in tqdm(sami_pass):
        sami_pass_min.append(i.dt.minute.values + 60*i.dt.hour.values + i.dt.second.values/60)
    interpolated_edens = []

    for i in tqdm(sat_final['sat_time']):
        sami_iter = []
        
        for j in sat_final['sami_time']:
            sami_iter.append(sat_final.e_dens.sel(sami_time=j, sat_time=i, 
                                                method='nearest').values)
        cs = CubicSpline(sami_pass_min, sami_iter)
        interpolated_edens.append(cs((i.dt.minute.values) + 60*(i.dt.hour.values)))

    sat_final['edens_sami_int'] = ('sat_time', interpolated_edens)
    return sat_final

    
def magnetic_coords(sat_data):
    alt = 850
    mlat = []
    mlon = []
    
    # Introducing the MLAT, MLON, MLT columns in satlocddf
    for glat, glon, gtime in tqdm(zip(sat_data.glat, sat_data.glon, sat_data.sat_time), total=len(sat_data.glon)):
        glon = glon.values
        glat = glat.values
        gtime = gtime.values
        date_str = (np.datetime_as_string(gtime))[:19]
        format_string = "%Y-%m-%dT%H:%M:%S"
        dtime = datetime.strptime(date_str, format_string)

        apex = Apex(dtime)
        magfield = apex.convert(glat, glon, 'geo', 'apex')        
        mlatitude, mlongitude = magfield
        mlat.append(mlatitude)
        mlon.append(mlongitude)

    # Adding mlat and mlon data variables in dataset
    sat_data['mlat'] = ('sat_time', mlat)
    sat_data['mlon'] = ('sat_time', mlon)
    
    return sat_data

