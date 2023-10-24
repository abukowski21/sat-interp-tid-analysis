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
from multiprocessing import Pool
from itertools import repeat


def interpolations(sat_info):
    '''
    The function defined to calculate the interpolations using 
    multiprocessing. The argument is a tuple of a tuple where the
    element tuple is made up of sami_time and edens (for fixed 
    sat_time) and sat_time at which edens need to be calculated by
    Cubic spline interpolation.
    '''
    sami_time, edens, sat_time = sat_info
    cs = CubicSpline(sami_time, edens)
    interpd_edens = cs(sat_time)
    
    
    return float(interpd_edens)


def edens_sattime(sat_data, day=None):
    """
    The purpose of the function is to do time interpolations such that 
    from the multiple Edens plots obtained for different Sami times
    we can get a single Edens plots that evolves with sat_step
    """

    if day is not None:
        sat_n = sat_data.where(sat_data.sat_time.dt.day == day, drop = 'True')
        sat_n = sat_n.where(sat_n.sami_time.dt.day == day, drop = 'True')
    else:
        sat_n = sat_data
    
    interpd_edens = np.zeros(sat_n.sat_step.shape)
    
    args = ()
    
    for t in tqdm(sat_n.sat_step.values):
    
        single_sat_step = sat_n.isel(sat_step=t)
        args1 = single_sat_step.sami_time
        args2 = single_sat_step.edens
        args3 = single_sat_step.sat_time
        args = args + ((args1,args2,args3),)

    interpd_edens = np.zeros(sat_n.sat_step.shape)
    
    with Pool(30) as pool:
    
        interpd_edens = pool.map(interpolations, args)
    
    return interpd_edens


def geo_to_mag_coord(geo_coords):
    '''
    Calculating the magnetic coordinates using parallel processing.
    The argument provided is a tuple made up of geographical 
    latitude, longitude, altitude, and time (a string of 
    "%Y-%m-%dT%H:%M:%S" form). 
    '''

    glat, glon, alt, gtime = geo_coords

    dtime_str = (np.datetime_as_string(gtime))[:19]
    format_string = "%Y-%m-%dT%H:%M:%S"
    dtime = datetime.strptime(dtime_str, format_string)
    
    decimal_year = dtime.year + ((dtime - datetime(dtime.year, 1, 1)).days)/365
    apex = Apex(decimal_year)
    mlat, mlon = apex.convert(glat, glon, 'geo', 'apex', height = alt)        
    m_coord = (mlat, mlon)

    return m_coord


def magnetic_coords(sat_data):
    '''
    Calculating the magnetic latitude and longitude for the satellite
    positions and add them in the satellite xarray variable
    '''

    glat = sat_data.glat.values
    glon = sat_data.glon.values
    alt = sat_data.alt.values
    gtime = sat_data.sat_time.values
    '''
    args = zip(glat, glon, alt, gtime)

    with Pool(30) as pool:

        mag_coords = pool.map(geo_to_mag_coord, args)
    
    mlat, mlon = mag_coords
    '''
    n = len(glat)
    mlat = np.zeros(n)
    mlon = np.zeros(n)
    
    for i in tqdm(range(n)):
        geo_coord = [glat[i], glon[i], alt[i], gtime[i]]
        m_coord = geo_to_mag_coord(geo_coord)    
        mlat[i], mlon[i] = m_coord
    
    sat_data['mlat'] = ('sat_step', mlat)
    sat_data['mlon'] = ('sat_step', mlon)
    
    return sat_data


def plot_edens(sat, nature):
    '''
    Plotting the characteristic given by the 'nature' variable for 6 
    different satellite passes. sat variable is having the data of 6
    satellite passes.
    '''
    
    # Plotting the edens on single plot. Variable is a string defining the data variable being plotted against mlat
    if len(sat) != 6:
        print('Error - Number of satellites in sat variable should be 6.\n')
        return
    
    color = ['black', 'purple', 'blue', 'yellow', 'orange', 'red']
    label = ['quiet 00:20 UTC', 'quiet 02:02 UTC', 'quiet 03:44 UTC', 'storm 15:39 UTC', 'storm 17:19 UTC', 'storm 19:01 UTC']

    plt.figure(figsize=(12,8))
    
    for s, c, l in zip(sat, color, label):
        plt.plot(s['glat'], s[nature], color=c, label=l)
        
    plt.legend(title='Satellite passes', loc = 'upper left')
    plt.xlabel('GLAT')
    plt.ylabel(nature)
    plt.title('Comparison Plots for Quiet and Storm Time Satellite Passes')
    plt.show()

    return
