import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import cartopy.crs as ccrs
from tqdm import tqdm
from apexpy import Apex
import os
import datetime as dt
import importlib as il

import sys
sys.path.append('../')
sys.path.append('/home/pxv220016/prasoon/data/sat_interp_repo/repo2/prasoon_utility_programs')
import functions


from p_tqdm import p_map
from multiprocessing import Pool
import itertools

month = 'sept'

if month == 'sept' or month == 'march':
    phase = 'equinox'
else:
    phase = 'solstice'

if month == 'sept':
    years = list(range(2010, 2024))
else:
    years = list(range(2021, 2025))


for year in tqdm(years):
    # Reading the data from output csv
    grnd_tec2 = pd.read_csv('/home/pxv220016/scratch/Qingyu_Cesar_EIA/outputs/' + month + '/' + str(year) + '_' + month + '_' + phase + '.csv')


    # Define the bin edges for MLAT and MLT
    bins_mlat = pd.cut(grnd_tec2['MLAT'], bins=pd.interval_range(start=-40, end=40, freq=1))
    bins_mlt = pd.cut(grnd_tec2['MLT'], bins=pd.interval_range(start=0, end=24, freq=0.25))
    
    # Create a new DataFrame with the bins
    grnd_tec2['MLAT_b'] = bins_mlat
    grnd_tec2['MLT_b'] = bins_mlt
    
    # Group by the bins (MLAT_b (primary) and MLT_b (secondary)) and calculate the average of TEC
    result = grnd_tec2.groupby(['MLT_b', 'MLAT_b'])['TEC'].mean().reset_index()
    # Converting the midpoint values of bins to float and assigning average TEC at those points
    result['MLAT_b'] = result['MLAT_b'].apply(lambda x: x.mid)
    result['MLT_b'] = result['MLT_b'].apply(lambda x: x.mid)
    result['MLAT_b'] = result['MLAT_b'].astype(float)
    result['MLT_b'] = result['MLT_b'].astype(float)


    result = result[result.MLT_b >= 5].reset_index(drop=True)
    
    
    filtered = []
    for t in result['MLT_b'].unique():
        result_f = result[result['MLT_b'] == t].reset_index(drop=True)
        fit = savgol_filter(np.array(result_f.TEC), 10, 2)
        filtered.extend(fit)
    result['TEC'] = filtered
    result = result.groupby(['MLAT_b', 'MLT_b']).sum().reset_index()
    
    
    # Identifying NH and SH peaks
    result_t = result[result.MLT_b >= 13].reset_index(drop=True)
    result_n = result_t[result_t.MLAT_b > 0].reset_index(drop=True)
    result_s = result_t[result_t.MLAT_b < 0].reset_index(drop=True)
    result_n = result_n.loc[result_n.groupby('MLT_b')['TEC'].idxmax()].reset_index(drop=True)
    result_s = result_s.loc[result_s.groupby('MLT_b')['TEC'].idxmax()].reset_index(drop=True)
    # Dropping the cases where the SH peak is not prominent and maxima appears to come at equator
    result_n.loc[result_n['MLAT_b'] < 8, 'TEC'] = np.nan
    result_s.loc[result_s['MLAT_b'] > -8, 'TEC'] = np.nan
    #print(result_s)
    
    i = 100*2*(result_n.TEC - result_s.TEC)/(result_n.TEC + result_s.TEC)
    result_ind = pd.DataFrame({'mlat_n': result_n.MLAT_b, 'tec_n': result_n.TEC, 'mlat_s': result_s.MLAT_b, 'tec_s': result_s.TEC, 'mlt': result_n.MLT_b, 'asy': i})
    result_ind.to_csv('/home/pxv220016/prasoon/data/sat_interp_repo/repo2/Qingyu_Cesar_EIA_IHA/outputs/' + month + '/asy_' + str(year) + '_' + month + '.csv', index=False)

    X, Y = np.meshgrid(result.MLT_b.unique(), result.MLAT_b.unique())
    Z = result.TEC.values.reshape(X.shape)
    contour_levels = list(range(0, 101, 10))

    # Creating new dataframe to remove rows where TEC= NaN
    result_n = pd.DataFrame({'mlt': result_ind.mlt, 'mlat_n': result_ind.mlat_n, 'tec_n': result_ind.tec_n})
    result_n = result_n.dropna(subset='tec_n')
    result_s = pd.DataFrame({'mlt': result_ind.mlt, 'mlat_s': result_ind.mlat_s, 'tec_s': result_ind.tec_s})
    result_s = result_s.dropna(subset='tec_s')

    
    fig = plt.figure(figsize=(12,4))
    specs = fig.add_gridspec(1, 2, width_ratios=[1,0.1])
    ax = []
    ax.append(fig.add_subplot(specs[0, 0]))
    c = ax[0].contourf(X, Y, Z, levels=contour_levels, cmap = 'jet')
    ax[0].scatter(result_n.mlt, result_n.mlat_n, c='b', s=30)
    ax[0].scatter(result_s.mlt, result_s.mlat_s, c='b', s=30)
    ax[0].set_title('Tracing EIA Peaks - ' + str(year) + ' ' + month + ' Equinox (Bins -> MLAT - 1 deg, MLT = 0.25 hour)')
    ax[0].set_ylabel('MLAT')
    ax[0].set_xlabel('MLT')
    ax[0].grid(True)
    
    cbar_ax = fig.add_subplot(specs[0,1])
    cbar = fig.colorbar(c, cax=cbar_ax, label='VTEC (TECUnits)', extend='both')
    fig.savefig('/home/pxv220016/prasoon/data/sat_interp_repo/repo2/Qingyu_Cesar_EIA_IHA/outputs/eia_peaks/'  + month + '/eia_peaks_'  + month + '_' + str(year) + '.jpg')
    fig.show()


