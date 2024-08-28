import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import math
import importlib as il
import datetime as dt
from scipy.interpolate import CubicSpline
from scipy.signal import find_peaks, savgol_filter


np.seterr(divide='ignore', invalid='ignore')


def ds_trim(ds, a, b, l1, l2):
    ds_n = ds.isel(Datetime=slice(a,b))
    ds_n = ds_n.where((ds_n.Glat > l1), drop=True)
    ds_n = ds_n.where((ds_n.Glat < l2), drop=True)
    return ds_n



def smoothen(func, lats):
    y = savgol_filter(func.dTEC.values, 5, 1)
    cs = CubicSpline(func.Glat.values, y)
    y1 = cs(lats)
    df = pd.DataFrame()
    df['Glat'] = lats
    df['dTEC'] = y1
    return df


def correlation(df1, df2):
    corr1 = []
    for i in range(len(df2)-10):
        l1 = len(df2)-10
        n = 0
        d1 = 0
        d2 = 0
        u1 = np.sum(df1[i:l1])/(l1-i)
        u2 = np.sum(df2[:l1-i])/(l1-i)
        
        for j,k in zip(df1[i:l1], df2[:l1-i]):
            n += (j-u1)*(k-u2)
            d1 += (j-u1)**2
            d2 += (k-u2)**2
        c = n/(np.sqrt(d1*d2))
        corr1.append(c)
    return corr1



def main(ds_n, l1, l2, date, show_speed = True, plot_3 = False, time1 = 0, time2 = 0, time3 = 0):
    speed = []
    lat_shift = []
    time = []
    
    i = 0
    j = 250
    
    while j < len(ds_n.dTEC.values):

        d = ds_n.Datetime.values[i]
        e = ds_n.Datetime.values[j]        
        
        df = ds_n.to_pandas()
        df_1 = (df.iloc[np.where((ds_n.Datetime == d))].reset_index(drop=True)).sort_values('Glat', ignore_index=True)
        df_2 = (df.iloc[np.where((ds_n.Datetime == e))].reset_index(drop=True)).sort_values('Glat', ignore_index=True)

        if len(df_1) < 10 or len(df_2) < 10:
            # Ensures that times with less data points are not used
            i += 250
            j += 250
            continue
        
        lats = np.linspace(l1,l2, 500)
        df1 = smoothen(df_1, lats)
        df2 = smoothen(df_2, lats)
        
        corr1 = correlation(df1.dTEC, df2.dTEC)
        m1 = max(corr1[:100])
        s1 = df2.Glat[corr1.index(m1)] - df2.Glat[0]
        t1 = int(e-d)/(60*10**9)
        slope1 = s1/t1
        
        speed.append(slope1)
        lat_shift.append(s1)
        time.append(ds_n.Datetime.values[int((i+j)/2)])
        i += 250
        j += 250

    ylim_max = int(max(speed) + 0.5)
    ylim_min = int(min(speed)) - 0.25
    speed_nonzero = [x for x in speed if x != 0]
    speed_av = sum(speed_nonzero)/len(speed_nonzero) # ignoring the zero terms of speed list
    
    # Plotting on Keogram
    fig1 = plt.figure(figsize=(12,8))
    specs1 = fig1.add_gridspec(1, 2, width_ratios=[1,0.1])
    ax1 = []
    ax1.append(fig1.add_subplot(specs1[0, 0]))
    c1 = ax1[0].scatter(ds_n.Datetime.values, ds_n.Glat.values, c=ds_n.dTEC.values, vmin=-1, vmax=1, cmap='viridis')
    ax1[0].set_title('dTEC Perturbations on GLAT-Time Plot (' + date + ')')
    ax1[0].set_xlabel('Time')
    ax1[0].set_ylabel('GLAT')   
    
    if show_speed:
        ax2 = ax1[0].twinx()
        ax2.plot(time, speed, color = 'r', linestyle = '--')
        ax2.scatter(time, speed, color = 'r')
        #ax2.set_ylim(ylim_min, ylim_max)
        ax2.axhline(y=speed_av, color='red', linestyle='--', label='Average Speed')
        ax2.set_ylabel('Speed (glat/min)', color = 'red')
        ax2.tick_params(axis='y', colors='red')
    
    if(plot_3):
        
        # Plotting 3 required lines
        ax1[0].axvline(x=ds_n.Datetime.values[time1], color='blue')
        ax1[0].axvline(x=ds_n.Datetime.values[time2], color='orange')
        ax1[0].axvline(x=ds_n.Datetime.values[time3], color='green')
        
        
        d = ds_n.Datetime.values[time1]
        e = ds_n.Datetime.values[time2]
        f = ds_n.Datetime.values[time3]
        
        df = ds_n.to_pandas()
        df_1 = (df.iloc[np.where((ds_n.Datetime == d))].reset_index(drop=True)).sort_values('Glat', ignore_index=True)
        df_2 = (df.iloc[np.where((ds_n.Datetime == e))].reset_index(drop=True)).sort_values('Glat', ignore_index=True)
        df_3 = (df.iloc[np.where((ds_n.Datetime == f))].reset_index(drop=True)).sort_values('Glat', ignore_index=True)
        
        df1 = smoothen(df_1, lats)
        df2 = smoothen(df_2, lats)
        df3 = smoothen(df_3, lats)
        
        plt.figure(figsize=(12,8))
        plt.plot(df1.Glat.values, df1.dTEC.values, label = str(d)[:-10])
        plt.plot(df2.Glat.values, df2.dTEC.values+2, label = str(e)[:-10])
        plt.plot(df3.Glat.values, df3.dTEC.values+4, label = str(f)[:-10])
        plt.title('Smooth (Savgol Ord=2) DTEC Values at Intervals of 10 Min')
        plt.xlabel('Glat')
        plt.ylabel('DTEC')
        plt.legend()


        corr1 = correlation(df1.dTEC, df2.dTEC)
        corr2 = correlation(df2.dTEC, df3.dTEC)
        
        m1 = max(corr1[:100])
        m2 = max(corr2[:100])
        s1 = df2.Glat[corr1.index(m1)] - df2.Glat[0]
        s2 = df3.Glat[corr2.index(m2)] - df3.Glat[0]
        print('Lat diff of peaks ' + str(s2-s1))
        
        print('Shift 1 =', s1, 'glats, Shift 2 =', s2, 'glats')
        t1 = int(e-d)/(60*10**9)
        slope1 = s1/t1
        
        t2 = int(f-e)/(60*10**9)
        slope2 = s2/t2
        print('slope 1 =', slope1, '(glat/min), Slope 2 =', slope2, '(glat/min)')

        
        plt.figure(figsize=(12,8))
        plt.plot(corr1, label = 'Corr_1 b/w ' + str(d)[11:-10] + ' and ' + str(e)[11:-10], color = 'yellow')
        plt.plot(corr2, label = 'Corr_2 b/w ' + str(e)[11:-10] + ' and ' + str(f)[11:-10], color = 'purple')
        plt.ylabel('Correlation')
        plt.xlabel('Shift (index)')
        plt.axvline(x=corr1.index(m1), color='yellow', linestyle='--', label='Max Corr_1')
        plt.axvline(x=corr2.index(m2), color='purple', linestyle='--', label='Max Corr_2')
        plt.title('Correlation for Plots for (' + str(d)[11:-10] + ' , ' + str(e)[11:-10] + ') and (' + str(e)[11:-10] + ' , ' + str(f)[11:-10] + ')')
        plt.legend()
        
        plt.figure(figsize=(12,8))
        plt.plot(df1.Glat.values, df1.dTEC.values, label = str(d)[:-10])
        plt.plot(df2.Glat.values + s1, df2.dTEC.values+2, label = str(e)[:-10])
        plt.plot(df3.Glat.values + s2 + s1, df3.dTEC.values+4, label = str(f)[:-10])
        plt.title('Shifted DTEC Values at 10 min Intervals (blu-org speed = ' + str(round(slope1, 4)) + ' glat/min, org-grn speed = ' + str(round(slope2, 4)) + ' glat/min)')
        plt.xlabel('Glat')
        plt.ylabel('DTEC')
        plt.legend()

    
    cbar_ax1 = fig1.add_subplot(specs1[0,1])
    cbar1 = fig1.colorbar(c1, cax=cbar_ax1, label='Perturbation TEC (TECUnits)', extend='both')
    fig1.show()

    return ds_n, speed, lat_shift