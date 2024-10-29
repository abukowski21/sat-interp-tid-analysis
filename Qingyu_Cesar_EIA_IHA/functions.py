import pandas as pd
import numpy as np
import datetime as dt


# FUNCTIONS DEFINITION

def day_to_date(day, year):
    leap = (year % 4 == 0) & ((year % 100 != 0) | (year % 400 == 0))
    days_in_month = np.array([31, 28 + leap, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    month = 0
    while day - days_in_month[month] > 0:
        day = day - days_in_month[month]
        month += 1
        if month == 12:
            month = 0
    return month + 1, day


def kp_index_filtering(sat_data, kp):
    # Here sat_data and kp are both dataframes
    # Kp values corresponding to date in sat_data is being stored in column of sat_data
    merged_data = pd.merge_asof(sat_data.sort_values('DT'), kp.sort_values('date'), left_on='DT', right_on='date', direction='backward')

    # Filter rows where 'kp' is less than 3
    filtered_data = merged_data[merged_data['kp'] <= 3].reset_index(drop=True)
    filtered_data = filtered_data.drop(['date'], axis = 1)
    return filtered_data


def magnetic_coords_parallel(sat_date, sat_glat, sat_glon, sat_tec):

    # Putting it here as apexpy can't be installed in conda env
    from apexpy import Apex
    
    gtime = sat_date
    decimal_year = gtime.year + ((gtime - dt.datetime(gtime.year, 1, 1)).days) / 365.25
    apex = Apex(decimal_year)
    
    mlat, mlon = apex.convert(sat_glat, sat_glon, 'geo', 'apex')
    mlt = apex.mlon2mlt(mlon, gtime)
    
    return [sat_date, sat_glat, sat_glon, sat_tec, mlat, mlon, mlt]


