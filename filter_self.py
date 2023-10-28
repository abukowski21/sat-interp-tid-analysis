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
import pyIGRF as pyigrf
import apexpy
from apexpy import Apex

from scipy.signal import butter, sosfiltfilt, sosfilt_zi, sosfilt, lfilter, filtfilt, savgol_filter, wiener, order_filter
from numpy.polynomial import Chebyshev
from moepy import lowess
from statsmodels.nonparametric.kernel_regression import KernelReg
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import fft
from scipy.fft import fftfreq


def FFT(y, sampling_freq=1):
    
    fourier = 2*fft(y)/len(y) # 2/len(y) is for normalization
    freq_axis = fftfreq(len(y), d=1.0/sampling_freq) # No need for Nyquist rate here
    
    return fourier, freq_axis


def filt_filt(da,
              freq=1/5,
              lims=[45, 80],
              order=1,
              xarray=True
    ):
    """
    Passing the Signal through a Band Pass Filter (BPF)

    Parameters
    ----------

    da : array_like or xarray dataset
        Data that needs to be filtered
    freq : float, optional
        Frequency of data collection (inverse of time taken to collect 
        two consecutive data). Default is 1/5.

        Example - for 1-second satellite data, `freq` = 1/1 sec^-1,
        similarly, for 5-minute satellite data, `freq` = 1/5 min^-1.
    lims : array_like, optional
        Critical time periods (or inverse of critical frequencies) 
        that contains the region that BPF will allow. The regions
        outside the lims will undergo perturbations. `lims` is a 
        length-2 sequence and its units should be similar to `freq`.

        Example - to filter the signals of period 150 to 300 seconds
        out of 1-second satellite data, lims = [150, 300], similarly,
        to filter the signal of period 45 to 80 minutes out of 
        5-minute satellite data, lims = [45, 80].
    order : int, optional
        Order of the Butterworth filter. Default is 1.
    xarray : bool, optional
        False if `da` is array_like , True is `da` is xarray dataset.
        Default is True.

    Returns
    -------
    filtd : np.array
        Filtered data
    filtd_perc : np.array
        Element-wise percentage difference from `da` data
    """

    if xarray:
        data = np.array(da.values)
    else:
        data = da

    # Defining critical frequencies as Wn
    lower_limit, upper_limit = lims
    Wn = np.array([1/upper_limit, 1/lower_limit])
    
    # Defining Nyquist rate
    nyquist_freq = 0.5*freq

    # Normalizing critical frequencies with Nyquist Rate: 
    Wn = Wn/nyquist_freq
    print(Wn)

    # Applying the Butterworth function to get Numerator (a) and
    # denominator (b) polynomials of IIR filter
    b, a = butter(order, Wn, btype='band', analog=False)

    # Passing b, a to filtfilt function to get filtered data
    filtd = np.array(filtfilt(b, a, data, axis=0))
    filtd_perc = np.array(100*(filtd)/data)

    return filtd, filtd_perc

'''
def filt_filt(da,
              sampling_freq=1,
              lims= [150, 300], #[0.0075, 0.02], #[0.0033, 0.0066],
              order=1,
              xarray=True):
    
    # Design the bandpass filter
    # lower_cutoff = 2*lower_lim/sampling_freq
    lower, upper = lims
    #nyquist_rate = 0.5*sampling_freq
        
    b, a = butter(order, [1/upper, 1/lower],
                  btype='band', analog=False, fs=sampling_freq)
    
    # Apply the filter to the data
    filtd = filtfilt(b, a, da, axis=0)

    if xarray:
        return filtd, (100*(filtd)/da.values)
    else:
        return filtd, (100*(filtd)/da)
'''

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

    
    