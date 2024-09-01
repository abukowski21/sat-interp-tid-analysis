import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import cartopy.crs as ccrs
from tqdm import tqdm
from apexpy import Apex
import os
import datetime as dt

import sys
sys.path.append('../')
sys.path.append('/home/pxv220016/prasoon/data/sat_interp_repo/repo2/prasoon_utility_programs')



def day_to_date(day, year):
    leap = (year % 4 == 0) & ((year % 100 != 0) | (year % 400 == 0))
    days_in_month = np.array([31, 28 + leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    month = 0
    while day - days_in_month[month] > 0:
        day = day - days_in_month[month]
        month += 1
    return month + 1, day


def kp_index_filtering(sat_data, kp):
    # Here sat_data and kp are both dataframes
    # Kp values corresponding to date in sat_data is being stored in column of sat_data
    merged_data = pd.merge_asof(sat_data.sort_values('DT'), kp, left_on='DT', right_on='date', direction='backward')

    # Filter rows where 'kp' is less than 3
    filtered_data = merged_data[merged_data['kp'] <= 3].reset_index(drop=True)
    filtered_data = filtered_data.drop(['date'], axis = 1)
    return filtered_data



def magnetic_coords(sat_data, apex):
    gtime = np.array(sat_data.DT)
    mlat, mlon = zip(*[apex.convert(i, j, 'geo', 'apex') for i, j in tqdm(zip(sat_data['GDLAT'], sat_data['GLON']))])
    mlt = apex.mlon2mlt(np.array(mlon), gtime)
    sat_data['MLAT'] = list(mlat)
    sat_data['MLON'] = list(mlon)
    sat_data['MLT'] = mlt
    
    return sat_data


def magnetic_coords_parallel(sat_date, sat_glat, sat_glon, apex):
    gtime = sat_date
    mlat, mlon = [apex.convert(sat_glat, sat_glon, 'geo', 'apex')]
    mlt = apex.mlon2mlt(mlon, gtime)
    
    return [sat_date, sat_glat, mlat, mlon, mlt]


