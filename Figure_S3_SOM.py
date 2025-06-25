#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 15:34:45 2025

@author: 6008399

Figure S3 SOM paper

"""

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import ruptures as rpt
from scipy.interpolate import CubicSpline
from scipy.interpolate import CubicHermiteSpline
import statsmodels.api as sm
import pandas as pd
from pandas.plotting import autocorrelation_plot
from pandas import DataFrame
from sklearn.linear_model import LinearRegression
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import xesmf as xe
import matplotlib.colors as mcolors
import cartopy.mpl.ticker as cticker
from scipy.optimize import fsolve
from numpy.polynomial.polynomial import Polynomial
from scipy.io import loadmat
from statsmodels.graphics.gofplots import qqplot

#Making pathway to folder with all data
directory = r'/../../Data/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%%

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_depth_0-1000m.nc','r')

time               = fh.variables['time'][:]         #time [model years]
som_star           = fh.variables['SOM'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_63-114_window_5_volume_integrated_SO30_finalcode.nc','r')
#fh = netcdf.Dataset(directory + 'TKE_year_63-114_window_5_volume_integrated_SO30.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_1_SO30 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_1_SO30 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_1_SO30 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_324-378_window_5_volume_integrated_SO30_finalcode.nc','r')
#fh = netcdf.Dataset(directory + 'TKE_year_324-378_window_5_volume_integrated_SO30.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 2
TKE_2_SO30 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_2_SO30 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_2_SO30 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_500-600_window_5_volume_integrated_SO30_finalcode.nc','r')
#fh = netcdf.Dataset(directory + 'TKE_year_500-600_window_5_volume_integrated_SO30.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 3
TKE_3_SO30 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_3_SO30 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_3_SO30 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

#fh = netcdf.Dataset(directory + 'KE_year_500-600_window_5_volume_integrated_SO30_fastcode.nc','r')
fh = netcdf.Dataset(directory + 'TKE_year_500-600_window_5_volume_integrated_SO30.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 3
TKE_3_SO30_Rene = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]

fh.close()

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data

#%%

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	lon		= fh.variables['lon'][:]
	lat		= fh.variables['lat'][:]
	time    = fh.variables['time'][:]
	TKE_int = fh.variables['TKE_int'][:]
	TKE		= fh.variables['TKE'][:]
	MKE		= fh.variables['MKE'][:]
	EKE		= fh.variables['EKE'][:]
	
	return lon, lat, time, TKE_int, TKE, MKE, EKE
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------


#Note that this TKE_int is not integrated over the whole Southern Ocean, but between -60 and 60E or so.
lon, lat, time_1, TKE_int_1, TKE_1, MKE_1, EKE_1	= ReadinData(directory+'KE_year_63-114_window_5_depth_averaged_5-317m_SO30.nc')
lon, lat, time_2, TKE_int_2, TKE_2, MKE_2, EKE_2	= ReadinData(directory+'KE_year_324-378_window_5_depth_averaged_5-317m_SO30.nc')
lon, lat, time_3, TKE_int_3, TKE_3, MKE_3, EKE_3	= ReadinData(directory+'KE_year_500-600_window_5_depth_averaged_5-317m_SO30.nc')

#%% Compare max and min TKE with each other

mean_TKE_1 = np.mean(TKE_int_1)
std_TKE_1  = np.std(TKE_int_1)

#Threshold for max phase
threshold_max = mean_TKE_1 + std_TKE_1
threshold_min = mean_TKE_1 - std_TKE_1

time_max  = time_1[TKE_int_1 > threshold_max]
EKE_1_max = EKE_1[TKE_int_1 > threshold_max]
TKE_1_max = TKE_1[TKE_int_1 > threshold_max]

time_min  = time_1[TKE_int_1 < threshold_min]
EKE_1_min = EKE_1[TKE_int_1 < threshold_min]
TKE_1_min = TKE_1[TKE_int_1 < threshold_min]

plt.figure()
plt.title('windows of max. and min. energetics we will look at SOM cycle 1')
plt.plot(time_1, TKE_1_SO30)
#plt.plot(EKE_3_SO30)

plt.axvline(x = time_max[0], color='red')
plt.axvline(x = time_max[-1], color='red')
plt.axvline(x = time_min[0], color='blue')
plt.axvline(x = time_min[-1], color='blue')

#%%Max and min second cycle

mean_TKE_2 = np.mean(TKE_int_2)
std_TKE_2  = np.std(TKE_int_2)

#Threshold for max phase
threshold_max = mean_TKE_2 + std_TKE_2
threshold_min = mean_TKE_2 - std_TKE_2

time_max  = time_2[TKE_int_2 > threshold_max]
EKE_2_max = EKE_2[TKE_int_2 > threshold_max]
TKE_2_max = TKE_2[TKE_int_2 > threshold_max]

time_min  = time_2[TKE_int_2 < threshold_min]
EKE_2_min = EKE_2[TKE_int_2 < threshold_min]
TKE_2_min = TKE_2[TKE_int_2 < threshold_min]

plt.figure()
plt.title('windows of max. and min. energetics we will look at SOM cycle 2')
plt.plot(time_2, TKE_2_SO30)
#plt.plot(EKE_3_SO30)

plt.axvline(x = time_max[0], color='red')
plt.axvline(x = time_max[-1], color='red')
plt.axvline(x = time_min[0], color='blue')
plt.axvline(x = time_min[-1], color='blue')

#%%Max and min third cycle

time_max  = time_3[32:53]
EKE_3_max = EKE_3[32:53]
TKE_3_max = TKE_3[32:53]

time_min  = time_3[67:88]
EKE_3_min = EKE_3[67:88]
TKE_3_min = TKE_3[67:88]

plt.figure()
plt.title('windows of max. and min. energetics we will look at SOM cycle 3')
plt.plot(time_3, TKE_3_SO30)
#plt.plot(EKE_3_SO30)

plt.axvline(x = time_max[0], color='red')
plt.axvline(x = time_max[-1], color='red')
plt.axvline(x = time_min[0], color='blue')
plt.axvline(x = time_min[-1], color='blue')

#%% EKE min/max subplots

# Create a figure with 2x2 subplots
fig, axes = plt.subplots(2, 2, subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(12, 10))

# Subplot 1: EKE difference regime A
ax1 = axes[0, 0]
CS1 = ax1.contourf(
    lon, lat, np.mean(EKE_2_min, axis=0) - np.mean(EKE_1_min, axis=0),
    levels=np.linspace(-15, 15.01, 31), extend='both', cmap='seismic', transform=ccrs.PlateCarree()
)
divider1 = make_axes_locatable(ax1)
ax_cb1 = divider1.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb1)
gl = ax1.gridlines(draw_labels=True)
gl.top_labels = False  
gl.right_labels = False  
gl.bottom_labels = True  
gl.left_labels = True  
cbar1 = plt.colorbar(CS1, cax=ax_cb1)
cbar1.set_label('$K_e$ difference [J]')
ax1.gridlines()
ax1.coastlines('50m')
ax1.add_feature(cfeature.LAND)
ax1.set_title('a) $K_e$ minimum difference (SOM 2 vs. 1)')

# Subplot 2: EKE difference regime B
ax2 = axes[0, 1]
CS2 = ax2.contourf(
    lon, lat, np.mean(EKE_2_max, axis=0) - np.mean(EKE_1_max, axis=0),
    levels=np.linspace(-15, 15.01, 31), extend='both', cmap='seismic', transform=ccrs.PlateCarree()
)
divider2 = make_axes_locatable(ax2)
ax_cb2 = divider2.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb2)
cbar2 = plt.colorbar(CS2, cax=ax_cb2)
cbar2.set_label('$K_e$ difference [J]')
gl = ax2.gridlines(draw_labels=True)
gl.top_labels = False  
gl.right_labels = False  
gl.bottom_labels = True  
gl.left_labels = True  
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND)
ax2.set_title('b) $K_e$ maximum difference  (SOM 2 vs. 1)')

# Subplot 3: Another variable (example)
ax3 = axes[1, 0]
CS3 = ax3.contourf(
    lon, lat, np.mean(EKE_3_min, axis=0) - np.mean(EKE_1_min, axis=0),
    levels=np.linspace(-15, 15.01, 31), extend='both', cmap='seismic', transform=ccrs.PlateCarree()
)
divider3 = make_axes_locatable(ax3)
ax_cb3 = divider3.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb3)
cbar3 = plt.colorbar(CS3, cax=ax_cb3)
cbar3.set_label('$K_e$ difference [J]')
gl = ax3.gridlines(draw_labels=True)
gl.top_labels = False  
gl.right_labels = False  
gl.bottom_labels = True  
gl.left_labels = True  
ax3.coastlines('50m')
ax3.add_feature(cfeature.LAND)
ax3.set_title('c) $K_e$ minimum difference  (SOM 3 vs. 1)')

# Subplot 4: Another variable (example)
ax4 = axes[1, 1]
CS4 = ax4.contourf(
    lon, lat, np.mean(EKE_3_max, axis=0) - np.mean(EKE_1_max, axis=0),
    levels=np.linspace(-15, 15.01, 31), extend='both', cmap='seismic', transform=ccrs.PlateCarree()
)
divider4 = make_axes_locatable(ax4)
ax_cb4 = divider4.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb4)
cbar4 = plt.colorbar(CS4, cax=ax_cb4)
cbar4.set_label('$K_e$ difference [J]')
gl = ax4.gridlines(draw_labels=True)
gl.top_labels = False  
gl.right_labels = False  
gl.bottom_labels = True  
gl.left_labels = True  
ax4.coastlines('50m')
ax4.add_feature(cfeature.LAND)
ax4.set_title('d) $K_e$ maximum difference  (SOM 3 vs. 1)')

# Adjust layout
plt.tight_layout()
plt.savefig(directory_figures + 'EKE_difference_min_max_3_2_1_0_318m.pdf')
plt.show()
