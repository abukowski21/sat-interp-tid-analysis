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
import pytz
from timezonefinder import TimezoneFinder


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


# Convert UTC to Local Time (Aaron's function)
def ut_to_lt(time_array, glon):
    """
    Compute local time from date and longitude.

    Parameters
    ----------
    time_array : array-like
        Array-like of datetime objects in universal time
    glon : array-like or float
        Float or array-like of floats containing geographic longitude
        in degrees. If single value or array of a different shape, all
        longitudes are applied to all times. If the shape is the same as
        `time_array`, the values are paired in the SLT calculation.

    Returns
    -------
    array of floats
        List of local times in hours
    """

    time_array = np.asarray(time_array)
    glon = np.asarray(glon)

    # Get UT seconds of day
    try:  # if numpy timestamps
        utsec = [(ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                  + ut.microsecond * 1.0e-6) / 3600.0 for ut in time_array]
    except BaseException:
        utsec = []
        for ut in time_array:
            ut = pd.Timestamp(ut)
            utsec.append((ut.hour * 3600.0 + ut.minute * 60.0 + ut.second
                          + ut.microsecond * 1.0e-6) / 3600.0)
    # Determine if the calculation is paired or broadcasted
    if glon.shape == time_array.shape:
        lt = np.array([utime + glon[i] / 15.0 for i,
                      utime in enumerate(utsec)])
    else:
        lt = np.array([utime + glon / 15.0 for utime in utsec])

    # Adjust to ensure that 0.0 <= lt < 24.0
    while np.any(lt < 0.0):
        lt[lt < 0.0] += 24.0

    while np.any(lt >= 24.0):
        lt[lt >= 24.0] -= 24.0
    
    return lt

'''
# Filter using FFT and then IFFT

def test_filt_filt(da,
              freq=1/5,
              lims=[45, 80],
              xarray=True,
              dimensions=1,
              size = 0
    ):
    print('yes')
    if xarray:
        data = np.array(da.values)
    else:
        data = da

    
    if dimensions == 2:
       
        print(len(data), len(data[0]), data[288][179], size, 'size comp')
        
        size2 = len(data[0])
        y = np.zeros([size, size2])
        y_perc = np.zeros([size, size2])
        
        num = 0
        
        for i in range(size2):
            d = np.zeros(size)
            for j in range(size):
                d[j] = data[j][i]
            fourier = np.fft.fft(d)
            freq_ax = np.fft.fftfreq(size, 1)

            Wn = [1.0/(freq*lims[1]), 1.0/(freq*lims[0])]
            fourier[np.where(np.logical_or(np.abs(freq_ax) < Wn[0], np.abs(freq_ax) > Wn[1]))] = 0
            y_iter = (np.fft.ifft(fourier)).real
            for j in range(size):
                y[j][i] = y_iter[j]
                y_perc[j][i] = 100*y[j][i]/d[j]
        return y, y_perc
        
    else:
        fourier = np.fft.fft(data)
        plt.plot(data)
        plt.show()
        freq_ax = np.fft.fftfreq(len(data), 1/freq)
        Wn = [1.0/(lims[1]), 1.0/(lims[0])]
        fourier[np.where(np.logical_or(np.abs(freq_ax) < Wn[0], np.abs(freq_ax) > Wn[1]))] = 0
        y = (np.fft.ifft(fourier)).real
        y_perc = [100*i/j for i, j in zip(y, data)]
        print(np.min(y_perc), np.max(y_perc))
    
        return y, y_perc
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
    yf = savgol_filter(y, 500, 3)
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

    
    