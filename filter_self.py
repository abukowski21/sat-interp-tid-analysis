import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import os, glob
import pandas as pd
from datetime import datetime, timedelta
from scipy.interpolate import CubicSpline
import importlib as il
from tqdm import tqdm

import pyIGRF as pyigrf
import apexpy
from apexpy import Apex

from scipy.signal import butter, sosfiltfilt, sosfilt_zi, sosfilt, lfilter, filtfilt, savgol_filter, wiener, order_filter
from numpy.polynomial import Chebyshev



def filt_filt(da,
              freq=5,
              lims=[40, 85],
              order=1,
              percent=True):

    '''
    custom filter function defined by Aaron as an example
    '''
    
    # Define sampling frequency and limits in minutes
    sampling_freq = freq
    lower_limit = min(lims)
    upper_limit = max(lims)

    # Convert limits to corresponding indices
    lower_index = int(lower_limit / sampling_freq)
    upper_index = int(upper_limit / sampling_freq)

    # Design the bandpass filter
    nyquist_freq = 0.5 * sampling_freq
    lower_cutoff = lower_index / nyquist_freq
    upper_cutoff = upper_index / nyquist_freq
    b, a = butter(order, [1/upper_cutoff, 1/lower_cutoff],
                  btype='band', analog=False)

    # Apply the filter to the data
    filtd = filtfilt(b, a, da, axis=0)
    # filtd = xr.apply_ufunc(filtfilt, b, a, da, dask='allowed')

    if percent:
        return (100*(filtd)/da)

    else:
        da.values = filtd
        return da



def polyfitting(degree, x, y):
    '''
    Filter defined for polynomial fitting. The arguments are degree 
    variable (degree of the polynomial you want to fit to), x (array of
    x data points) and y (input data points that need to be filtered or
    smoothened) 
    '''
    
    coeff = np.polyfit(x,y,degree)
    yf = np.zeros(len(x))
    m = 0
    for i in x:
        sum = 0
        k = 0
        for k in range(degree+1):            
            sum += coeff[k]*(pow(i, degree-k))
        yf[m] = sum
        m += 1
    return yf


def best_filters(x,y):
    yfit = []
    name = []

    '''
    # order filter - Not good
    yf = order_filter(y, np.ones(101), rank=50)
    yfit.append(yf)
    name.append('Order Filter')
    '''

    '''
    # wiener filter - not good
    yf = wiener(y,180)
    yfit.append(yf)
    name.append('Wiener Filter')
    '''

    # Savitzky-Golay filter
    yf = savgol_filter(y, 500, 4)
    yfit.append(yf)
    name.append('Savgol Filter')

    # Polyfitting
    yf = polyfitting(10, x, y)
    yfit.append(yf)
    name.append('Polynomial Fitting')

    # Chebyshev Fitting
    c = Chebyshev.fit(x, y, 10)
    yf = c(x)
    yfit.append(yf)
    name.append('Chebyshev Fitting')

    # Rolling mean
    window = 50
    yf1 = np.convolve(y, np.ones(window)/window, mode='valid')
    yf = np.zeros(len(y))
    for i in range(len(y)):
        if i-1 < window:
            yf[i] = np.NaN
        
        else:
            yf[i] = yf1[i-window]    

    yfit.append(yf)
    name.append('Rolling mean')
    
    return name, yfit



    
    