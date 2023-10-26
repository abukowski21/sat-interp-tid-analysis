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
from pandas import DataFrame
import math

from scipy.signal import butter, sosfiltfilt, sosfilt_zi, sosfilt, lfilter, filtfilt, savgol_filter, wiener, order_filter
from numpy.polynomial import Chebyshev
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft
from scipy.fft import fftfreq


def FFT(y, sampling_freq=1):
    
    fourier = 2*fft(y)/len(y) # 2/len(y) is for normalization
    freq_axis = fftfreq(len(y), d=1.0/sampling_freq) # No need for Nyquist rate here
    
    return fourier, freq_axis

def ab_filt_filt(da,
              freq=5,
              lims=[45, 80],
              order=1,
              xarray=True):

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
    filtd = np.array(filtfilt(b, a, da, axis=0))
    # filtd = xr.apply_ufunc(filtfilt, b, a, da, dask='allowed')

    if xarray:
        return filtd, (100*(filtd)/da.values)
    else:
        return filtd, (100*(filtd)/da)


def filt_filt(da,
              sampling_freq=1,
              lims= [150, 300], #[0.0075, 0.02], #[0.0033, 0.0066],
              order=1,
              xarray=True):
    
    # Design the bandpass filter
    # lower_cutoff = 2*lower_lim/sampling_freq
    lower, upper = lims
    nyquist_rate = 0.5*sampling_freq
    #print('lower limit =', 1/upper, ', upper limit =', 1/lower)
    
    # 1/lower will give the upper-frequency limit of bandpass and vice versa
    b, a = butter(order, [1/upper, 1/lower],
                  btype='band', analog=False, fs=nyquist_rate)
    
    # Apply the filter to the data
    filtd = filtfilt(b, a, da, axis=0)

    if xarray:
        return filtd, (100*(filtd)/da.values)
    else:
        return filtd, (100*(filtd)/da)

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

def rolling_mean(y, window):
    
    ydf = pd.DataFrame(y, columns=["edens"])
    yf = ydf['edens'].rolling(window = window, center=True).mean()
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
    yf = rolling_mean(y, 180)
    yfit.append(yf)
    name.append('Rolling Mean')

    '''
    # Custom filter
    yf = filt_filt(y)
    yfit.append(y-yf)  # Sending the perturbations CHANGE LATER
    name.append('Custom Filter')
    '''
    
    '''
    # Triangular Moving Average
    yf = rolling_mean(y, 100)
    yf = rolling_mean(yf, 50)
    yfit.append(yf)
    name.append('Triangular Moving Average')
    '''

    # LOWESS Filter
    from moepy import lowess
    lowess_m = lowess.Lowess()
    lowess_m.fit(x,y,0.17)
    yf = lowess_m.predict(x)
    yfit.append(yf)
    name.append('LOWESS Fit')

    
    return name, yfit


def DFT(x):
    """
    Function to calculate the 
    discrete Fourier Transform 
    of a 1D real-valued signal x
    """

    N = len(x)
    n = np.arange(N)
    k = n.reshape((N, 1))
    e = np.exp(-2j * np.pi * k * n / N)
    
    X = np.dot(e, x)
    
    return X

    
    