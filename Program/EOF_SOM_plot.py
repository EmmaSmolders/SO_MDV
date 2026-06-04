#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 21:22:20 2025

@author: 6008399
"""

#%%

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
from scipy.interpolate import CubicSpline
from scipy.interpolate import CubicHermiteSpline
import statsmodels.api as sm
import pandas as pd
from pandas.plotting import autocorrelation_plot
from pandas import DataFrame
import cartopy.crs as ccrs
import numpy as np
import xarray as xr
import matplotlib.colors as mcolors
import cartopy.mpl.ticker as cticker

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/HR_POP/Zenodo/Data_Final/'
directory_figures	= '/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%% Read in data

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	lon 		= fh.variables['lon'][:]			#Longitude
	lat 		= fh.variables['lat'][:]			#Latitude 
	eof         = fh.variables['eof'][:]           #number of EOFs
	time		= fh.variables['time'][:]			#Model year
	PC		    = fh.variables['PC'][:] 			#Principal component
	VAR		    = fh.variables['VAR'][:]	 		#Variance of the PCs/EOFs
	EOF	       	= fh.variables['EOF'][:]			#EOFs

	fh.close()

	return lon, lat, eof, time, PC, VAR, EOF

#%%

moving_average = 5
month_start    = 1
month_end      = 12

region = 'SO' #region WG, AU or SO, or Pacific

lon, lat, eof_E1, time_E1, PC_E1, VAR_E1, EOF_E1		= ReadinData(directory + 'EOF_SOM_SST_month_1_12_moving_average_5_POP_year_3_98.nc')
lon, lat, eof_E2, time_E2, PC_E2, VAR_E2, EOF_E2		= ReadinData(directory + 'EOF_SOM_SST_month_1_12_moving_average_5_POP_year_502_598.nc')

#%%

plt.figure()
plt.contourf(EOF_E1[0,:,:], levels = np.linspace(-0.02,0.02,21), cmap='RdBu_r')

#%% Take first EOF for NAO

EOF_TEMP_E1 = EOF_E1[0,:,:]
EOF_TEMP_E2 = EOF_E2[0,:,:]


#%%

#Align signs using correlation of PCs (first mode)
corr = np.corrcoef(EOF_TEMP_E1, EOF_TEMP_E2)[0,1]
if corr < 0:
    print('Switching signs')
    EOF_TEMP_E2 *= -1
    PC_E2[0,:] *= -1
    
#EOF_SLP_E4 = -EOF_SLP_E4
#PC_E4 = -PC_E4

#%%    

fig, axs = plt.subplots(2, 2, figsize=(11, 6), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})

for ax in axs.flat:  
    ax.coastlines()
    
plt.suptitle('SST '+str(region)+'', fontsize=15)
    
c1 = axs[0,0].contourf(lon, lat, EOF_TEMP_E1, transform=ccrs.PlateCarree(), levels = np.linspace(-0.02,0.02,21), cmap='RdBu_r')
fig.colorbar(c1, ax=axs[0,0], orientation='horizontal')
axs[0,0].set_title('a) First EOF AMOC on (var.ex. = '+str(int(VAR_E1[0]))+'%)')
axs[0,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,0].xaxis.set_major_formatter(lon_formatter)
axs[0,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,0].yaxis.set_major_formatter(lat_formatter)
#axs[0,0].set_xlim(-10, 120)
#axs[0,0].set_ylim(-70, -35)

c2 = axs[0,1].contourf(lon, lat, EOF_TEMP_E2, transform=ccrs.PlateCarree(), levels = np.linspace(-.02,0.02,21), cmap='RdBu_r')
fig.colorbar(c2, ax=axs[0,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[0,1].set_title('b) First EOF AMOC off (var.ex. = '+str(int(VAR_E2[0]))+'%)')
axs[0,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,1].xaxis.set_major_formatter(lon_formatter)
axs[0,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,1].yaxis.set_major_formatter(lat_formatter)
#axs[0,1].set_xlim(-10, 120)
#axs[0,1].set_ylim(-70, -35)

#Define region for som regions
lon_min, lon_max = 175, 220
lat_min, lat_max = -60, -45

x_1	= np.arange(lon_min, lon_max + 0.1, 0.1)
y_1	= np.zeros(len(x_1)) + lat_min
y_2	= np.zeros(len(x_1)) + lat_max

y_3 = np.arange(lat_min, lat_max + 0.1, 0.1)
x_2 = np.zeros(len(y_3)) + lon_min
x_3 = np.zeros(len(y_3)) + lon_max

axs[0,1].plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

axs[1,1].plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1,1].plot(x_1, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1,1].plot(x_2, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1,1].plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

c1 = axs[1,0].contourf(lon, lat, EOF_E1[1,:,:], transform=ccrs.PlateCarree(), levels = np.linspace(-0.02,0.02,21), cmap='RdBu_r')
fig.colorbar(c1, ax=axs[1,0], orientation='horizontal')
axs[1,0].set_title('c) Second EOF AMOC on (var.ex. = '+str(int(VAR_E1[1]))+'%)')
axs[1,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,0].xaxis.set_major_formatter(lon_formatter)
axs[1,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,0].yaxis.set_major_formatter(lat_formatter)
axs[1,0].set_xlim(-10, 120)
axs[1,0].set_ylim(-70, -35)

c2 = axs[1,1].contourf(lon, lat, EOF_E2[1,:,:], transform=ccrs.PlateCarree(), levels = np.linspace(-.02,0.02,21), cmap='RdBu_r')
fig.colorbar(c2, ax=axs[1,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[1,1].set_title('d) Second EOF AMOC off (var.ex. = '+str(int(VAR_E2[1]))+'%)')
axs[1,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,1].xaxis.set_major_formatter(lon_formatter)
axs[1,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,1].yaxis.set_major_formatter(lat_formatter)
axs[1,1].set_xlim(-10, 120)
axs[1,1].set_ylim(-70, -35)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'EOF_SST_'+str(region)+'_moving_average_'+str(moving_average)+'_HRPOP.pdf')
plt.show()    


#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=0)})

for ax in axs.flat:
    ax.coastlines()

c1 = axs[0,0].contourf(lon, lat, EOF_TEMP_E2 - EOF_TEMP_E1, transform=ccrs.PlateCarree(), levels = np.linspace(-0.02,0.02,21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[0,0], orientation='horizontal')
axs[0,0].set_title('a) Difference first EOF pattern')
axs[0,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,0].xaxis.set_major_formatter(lon_formatter)
axs[0,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,0].yaxis.set_major_formatter(lat_formatter)
#axs[0,0].set_xlim(-50, 50)
axs[0,0].set_ylim(-70, -35)

c2 = axs[0,1].contourf(lon, lat, EOF_E2[1,:,:] - EOF_E1[1,:,:], transform=ccrs.PlateCarree(), levels = np.linspace(-.02,0.02,21), cmap='RdBu_r', extend = 'both')
fig.colorbar(c2, ax=axs[0,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[0,1].set_title('b) Difference second EOF pattern')
axs[0,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,1].xaxis.set_major_formatter(lon_formatter)
axs[0,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,1].yaxis.set_major_formatter(lat_formatter)
#axs[0,1].set_xlim(-50, 50)
axs[0,1].set_ylim(-70, -35)

c1 = axs[1,0].contourf(lon, lat, EOF_TEMP_E2 * np.max(abs(PC_E2[0,:])) - EOF_TEMP_E1 * np.max(abs(PC_E1[0,:])), transform=ccrs.PlateCarree(), levels = np.linspace(-0.003,0.003,21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[1,0], orientation='horizontal')
axs[1,0].set_title('c) Difference first EOF*max(|PC|)')
axs[1,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,0].xaxis.set_major_formatter(lon_formatter)
axs[1,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,0].yaxis.set_major_formatter(lat_formatter)
#axs[1,0].set_xlim(-50, 50)
axs[1,0].set_ylim(-70, -35)

c2 = axs[1,1].contourf(lon, lat, EOF_E2[1,:,:] * np.max(abs(PC_E2[1,:])) - EOF_E1[1,:,:] * np.max(abs(PC_E1[1,:])), transform=ccrs.PlateCarree(), levels = np.linspace(-.003,0.003,21), cmap='RdBu_r', extend = 'both')
fig.colorbar(c2, ax=axs[1,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[1,1].set_title('d) Difference second EOF*max(|PC|)')
axs[1,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,1].xaxis.set_major_formatter(lon_formatter)
axs[1,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,1].yaxis.set_major_formatter(lat_formatter)
#axs[1,1].set_xlim(-50, 50)
axs[1,1].set_ylim(-70, -35)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'EOF_SST_'+str(region)+'_diff_moving_average_'+str(moving_average)+'_HRPOP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(11, 6), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})

for ax in axs.flat:  
    ax.coastlines()
    
plt.suptitle('SST '+str(region)+'', fontsize=15)
    
c1 = axs[0,0].contourf(lon, lat, EOF_TEMP_E1, transform=ccrs.PlateCarree(), levels = np.linspace(-0.02,0.02,21), cmap='RdBu_r')
fig.colorbar(c1, ax=axs[0,0], orientation='horizontal')
axs[0,0].set_title('a) First EOF AMOC on (var.ex. = '+str(int(VAR_E1[0]))+'%)')
axs[0,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,0].xaxis.set_major_formatter(lon_formatter)
axs[0,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,0].yaxis.set_major_formatter(lat_formatter)
axs[0,0].set_xlim(-10, 120)
axs[0,0].set_ylim(-70, -35)

c2 = axs[0,1].contourf(lon, lat, EOF_TEMP_E2, transform=ccrs.PlateCarree(), levels = np.linspace(-.02,0.02,21), cmap='RdBu_r')
fig.colorbar(c2, ax=axs[0,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[0,1].set_title('b) First EOF AMOC off (var.ex. = '+str(int(VAR_E2[0]))+'%)')
axs[0,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,1].xaxis.set_major_formatter(lon_formatter)
axs[0,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,1].yaxis.set_major_formatter(lat_formatter)
axs[0,1].set_xlim(-10, 120)
axs[0,1].set_ylim(-70, -35)

axs[0,1].plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

axs[1,1].plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1,1].plot(x_1, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1,1].plot(x_2, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1,1].plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

c1 = axs[1,0].contourf(lon, lat, EOF_TEMP_E2 * np.max(abs(PC_E2[0,:])) - EOF_TEMP_E1 * np.max(abs(PC_E1[0,:])), transform=ccrs.PlateCarree(), levels = np.linspace(-0.005,0.005,21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[1,0], orientation='horizontal')
axs[1,0].set_title('c) Difference first EOF*max(|PC|)')
axs[1,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,0].xaxis.set_major_formatter(lon_formatter)
axs[1,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,0].yaxis.set_major_formatter(lat_formatter)
axs[1,0].set_xlim(-10, 120)
axs[1,0].set_ylim(-70, -35)

c1 = axs[1,1].contourf(lon, lat, EOF_TEMP_E2 - EOF_TEMP_E1, transform=ccrs.PlateCarree(), levels = np.linspace(-0.04,0.04,21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[1,1], orientation='horizontal')
axs[1,1].set_title('d) Difference first EOF pattern')
axs[1,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,1].xaxis.set_major_formatter(lon_formatter)
axs[1,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,1].yaxis.set_major_formatter(lat_formatter)
axs[1,1].set_xlim(-10, 120)
axs[1,1].set_ylim(-70, -35)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'EOF1_diff_SST_'+str(region)+'_moving_average_'+str(moving_average)+'_HRPOP.pdf')
plt.show()    


#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=0)})

for ax in axs.flat:
    ax.coastlines()

c1 = axs[0,0].contourf(lon, lat, EOF_TEMP_E2 - EOF_TEMP_E1, transform=ccrs.PlateCarree(), levels = np.linspace(-0.02,0.02,21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[0,0], orientation='horizontal')
axs[0,0].set_title('a) Difference first EOF pattern')
axs[0,0].set_xticks(np.arange(-50,50,20), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,0].xaxis.set_major_formatter(lon_formatter)
axs[0,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,0].yaxis.set_major_formatter(lat_formatter)
axs[0,0].set_xlim(-80, 80)
axs[0,0].set_ylim(-70, -35)

c2 = axs[0,1].contourf(lon, lat, EOF_E2[1,:,:] - EOF_E1[1,:,:], transform=ccrs.PlateCarree(), levels = np.linspace(-.02,0.02,21), cmap='RdBu_r', extend = 'both')
fig.colorbar(c2, ax=axs[0,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[0,1].set_title('b) Difference second EOF pattern')
axs[0,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[0,1].xaxis.set_major_formatter(lon_formatter)
axs[0,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[0,1].yaxis.set_major_formatter(lat_formatter)
#axs[0,1].set_xlim(-50, 50)
axs[0,1].set_ylim(-70, -35)

c1 = axs[1,0].contourf(lon, lat, EOF_TEMP_E2 * np.max(abs(PC_E2[0,:])) - EOF_TEMP_E1 * np.max(abs(PC_E1[0,:])), transform=ccrs.PlateCarree(), levels = np.linspace(-0.003,0.003,21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[1,0], orientation='horizontal')
axs[1,0].set_title('c) Difference first EOF*max(|PC|)')
axs[1,0].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,0].xaxis.set_major_formatter(lon_formatter)
axs[1,0].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,0].yaxis.set_major_formatter(lat_formatter)
#axs[1,0].set_xlim(-50, 50)
axs[1,0].set_ylim(-70, -35)

c2 = axs[1,1].contourf(lon, lat, EOF_E2[1,:,:] * np.max(abs(PC_E2[1,:])) - EOF_E1[1,:,:] * np.max(abs(PC_E1[1,:])), transform=ccrs.PlateCarree(), levels = np.linspace(-.003,0.003,21), cmap='RdBu_r', extend = 'both')
fig.colorbar(c2, ax=axs[1,1], orientation='horizontal')
#axs[1].set_title('b) Second EOF (var.ex. = '+str(int(VAR[1]))+'%)')
axs[1,1].set_title('d) Difference second EOF*max(|PC|)')
axs[1,1].set_xticks(np.arange(-180,180,50), crs=ccrs.PlateCarree())
lon_formatter = cticker.LongitudeFormatter()
axs[1,1].xaxis.set_major_formatter(lon_formatter)
axs[1,1].set_yticks(np.arange(-80,-30,20), crs=ccrs.PlateCarree())
lat_formatter = cticker.LatitudeFormatter()
axs[1,1].yaxis.set_major_formatter(lat_formatter)
#axs[1,1].set_xlim(-50, 50)
axs[1,1].set_ylim(-70, -35)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'EOF_SST_'+str(region)+'_diff_moving_average_'+str(moving_average)+'_HRPOP.pdf')
plt.show()

#%% plot PC's

#Central moving average
def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

window = 20

fig, axs = plt.subplots(1, 2, figsize=(14, 4))  

plt.suptitle('First PC patterns', fontsize=14)

axs[0].set_title('a) Model year 1-100')
axs[0].plot(time_E1, PC_E1[0,:], color='red', alpha = 0.3, label='AMOC on')
axs[0].plot(time_E1[window//2 : -window//2 + 1], Moving_average(PC_E1[0,:], window), color='red')
axs[0].set_ylim(-0.2, 0.35)
axs[0].legend()

axs[1].set_title('b) Model year 500-600')
axs[1].plot(time_E2, PC_E2[0,:], color='red', alpha = 0.3, label='AMOC off')
axs[1].plot(time_E2[window//2 : -window//2 + 1], Moving_average(PC_E2[0,:], window), color='red')
axs[1].set_ylim(-0.2, 0.35)
axs[1].legend()

fig, axs = plt.subplots(1, 2, figsize=(14, 4))  
plt.suptitle('Second PC patterns', fontsize=14)

axs[0].set_title('a) Model year 1-100')
axs[0].plot(time_E1, PC_E1[1,:], color='blue', alpha = 0.3, label='AMOC on')
axs[0].plot(time_E1[window//2 : -window//2 + 1], Moving_average(PC_E1[1,:], window), color='blue')
axs[0].set_ylim(-0.2, 0.35)
axs[0].legend()

axs[1].set_title('b) Model year 500-600')
axs[1].plot(time_E2, PC_E2[1,:], color='blue', alpha = 0.3, label='AMOC off')
axs[1].plot(time_E2[window//2 : -window//2 + 1], Moving_average(PC_E2[1,:], window), color='blue')
axs[1].set_ylim(-0.2, 0.35)
axs[1].legend()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6)) 
# First PC patterns
axs[0, 0].set_title('a) First PC: AMOC on')
axs[0, 0].plot(time_E1, PC_E1[0, :], color='red', alpha=0.3, label='AMOC on')
axs[0, 0].plot(time_E1[window // 2: -window // 2 + 1], Moving_average(PC_E1[0, :], window), color='red')
axs[0, 0].set_ylim(-0.2, 0.35)
axs[0, 0].legend()

axs[0, 1].set_title('b) First PC: AMOC off')
axs[0, 1].plot(time_E2, PC_E2[0, :], color='red', alpha=0.3, label='AMOC off')
axs[0, 1].plot(time_E2[window // 2: -window // 2 + 1], Moving_average(PC_E2[0, :], window), color='red')
axs[0, 1].set_ylim(-0.2, 0.35)
axs[0, 1].legend()

# Second PC patterns
axs[1, 0].set_title('c) Second PC: AMOC on')
axs[1, 0].plot(time_E1, PC_E1[1, :], color='blue', alpha=0.3, label='AMOC on')
axs[1, 0].plot(time_E1[window // 2: -window // 2 + 1], Moving_average(PC_E1[1, :], window), color='blue')
axs[1, 0].set_ylim(-0.2, 0.35)
axs[1, 0].legend()

axs[1, 1].set_title('d) Second PC: AMOC off')
axs[1, 1].plot(time_E2, PC_E2[1, :], color='blue', alpha=0.3, label='AMOC off')
axs[1, 1].plot(time_E2[window // 2: -window // 2 + 1], Moving_average(PC_E2[1, :], window), color='blue')
axs[1, 1].set_ylim(-0.2, 0.35)
axs[1, 1].legend()

plt.tight_layout()
plt.savefig(directory_figures + 'PCs_AMOC_off_on_HRPOP_'+str(region)+'.pdf')
plt.show()

#%%

import cartopy.feature as cfeature



# Create a 1x2 subplot with South Polar Stereographic projection
fig = plt.figure(figsize=(10, 8))

axs = np.empty((2,2), dtype=object)

# Top row — South Polar maps
axs[0,0] = fig.add_subplot(2, 2, 1, projection=ccrs.SouthPolarStereo())
axs[0,1] = fig.add_subplot(2, 2, 2, projection=ccrs.SouthPolarStereo())

# Bottom row — normal (Cartesian) plots
axs[1,0] = fig.add_subplot(2, 2, 3)
axs[1,1] = fig.add_subplot(2, 2, 4)

# Add coastlines and gridlines to each subplot
axs[0,0].coastlines()
axs[0,0].set_extent([-180, 180, -90, -35], crs=ccrs.PlateCarree())  # Focus on Antarctica

axs[0,1].coastlines()
axs[0,1].set_extent([-180, 180, -90, -35], crs=ccrs.PlateCarree())  # Focus on Antarctica

    # Add gridlines
gl = axs[0,0].gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

gl = axs[0,1].gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

axs[0,0].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none')

axs[0,1].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none')

c1 = axs[0,0].contourf(lon, lat, -EOF_TEMP_E1, transform=ccrs.PlateCarree(), levels=np.linspace(-0.025, 0.025, 21), cmap='RdBu_r', extend='both')
fig.colorbar(c1, ax=axs[0,0], orientation='vertical', shrink = 1)
axs[0,0].set_title('a) First EOF - AMOC on (var.ex. = '+str(int(VAR_E1[0]))+'%)', fontsize=13)

# Second subplot: EOF_SLP_E4
c2 = axs[0,1].contourf(lon, lat, -EOF_TEMP_E2, transform=ccrs.PlateCarree(), levels=np.linspace(-0.025, 0.025, 21), cmap='RdBu_r', extend='both')
fig.colorbar(c2, ax=axs[0,1], orientation='vertical', shrink=1)
axs[0,1].set_title('b) First EOF - AMOC off (var.ex. = '+str(int(VAR_E2[0]))+'%)', fontsize=13)


#Define region for som regions
lon_min, lon_max = -50, 0
lat_min, lat_max = -50, -35

lon_min_star, lon_max_star = 150, 190
lat_min_star, lat_max_star = -70, -60

lon_min_p, lon_max_p = 170, 210
lat_min_p, lat_max_p = -60, -45

#ACC transect
acc_lon_min, acc_lon_max = -68, -64
acc_lat_min, acc_lat_max = -66.5, -55

#WGKP
lon_min_wgkp, lon_max_wgkp = -35, 80
lat_min_wgkp, lat_max_wgkp = -80, -50

x_1	= np.arange(lon_min, lon_max + 0.1, 0.1)
y_1	= np.zeros(len(x_1)) + lat_min
y_2	= np.zeros(len(x_1)) + lat_max

y_3 = np.arange(lat_min, lat_max + 0.1, 0.1)
x_2 = np.zeros(len(y_3)) + lon_min
x_3 = np.zeros(len(y_3)) + lon_max

axs[0,0].plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,0].plot(x_1, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,0].plot(x_2, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,0].plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_star	= np.arange(lon_min_star, lon_max_star + 0.1, 0.1)
y_1_star	= np.zeros(len(x_1_star)) + lat_min_star
y_2_star	= np.zeros(len(x_1_star)) + lat_max_star

y_3_star = np.arange(lat_min_star, lat_max_star + 0.1, 0.1)
x_2_star = np.zeros(len(y_3_star)) + lon_min_star
x_3_star = np.zeros(len(y_3_star)) + lon_max_star

x_1_p	= np.arange(lon_min_p, lon_max_p + 0.1, 0.1)
y_1_p	= np.zeros(len(x_1_p)) + lat_min_p
y_2_p	= np.zeros(len(x_1_p)) + lat_max_p

y_3_p = np.arange(lat_min_p, lat_max_p + 0.1, 0.1)
x_2_p = np.zeros(len(y_3_p)) + lon_min_p
x_3_p = np.zeros(len(y_3_p)) + lon_max_p

axs[0,1].plot(x_1_p, y_1_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_p, y_2_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_p, y_3_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_p, y_3_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)


x_1_acc	= np.arange(acc_lon_min, acc_lon_max + 0.1, 0.1)
y_1_acc	= np.zeros(len(x_1_acc)) + acc_lat_min
y_2_acc	= np.zeros(len(x_1_acc)) + acc_lat_max

y_3_acc = np.arange(acc_lat_min, acc_lat_max + 0.1, 0.1)
x_2_acc = np.zeros(len(y_3_acc)) + -66
x_3_acc = np.zeros(len(y_3_acc)) + acc_lon_max

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp

axs[1,0].set_title('c) First PC - AMOC on', fontsize=13)
axs[1,0].plot(time_E1, -PC_E1[0,:], color='darkgreen')
#axs[1,0].plot(time_E1[window//2 : -window//2 + 1], Moving_average(-PC_E1[0,:], window), color='darkgreen', label='Smoothed')
axs[1,0].set_ylim(-0.2, 0.2)
#axs[1,0].legend()
axs[1,0].grid()
axs[1,0].set_xlabel('Time [model year]', fontsize=12)


axs[1,1].set_title('d) First PC - AMOC off', fontsize=13)
axs[1,1].plot(time_E2, -PC_E2[0,:], color='darkgreen')
#axs[1,1].plot(time_E2[window//2 : -window//2 + 1], Moving_average(-PC_E2[0,:], window), color='darkgreen', label='Smoothed')
axs[1,1].set_ylim(-0.2, 0.2)
##axs[1,1].legend()
axs[1,1].grid()
axs[1,1].set_xlabel('Time [model year]', fontsize=12)


# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'Figure_2_SOM_OS.pdf')
plt.show()




# %% North et al. (1982) test for significance of EOFs

import numpy as np

def north_test(var, time, mode1=0, mode2=1):
    """
    North et al. (1982) rule-of-thumb test for separation of EOF modes.

    Parameters
    ----------
    var : array-like
        Eigenvalues or explained variance percentages of EOF modes. It's usually done with eigenvalues, but the explained variance are just the eigenvalues normalized by their sum, so it should work as well.
    time : array-like
        Time axis used in the EOF analysis.
    mode1, mode2 : int
        Indices of the modes to compare (Python indexing).

    Returns
    -------
    dict with test results
    """
    var = np.asarray(var)
    N = len(time)
    #N = effective_N(time)

    lam1 = var[mode1]
    lam2 = var[mode2]

    delta1 = lam1 * np.sqrt(2 / N)
    delta2 = lam2 * np.sqrt(2 / N)

    separated_simple = (lam1 - lam2) > delta1
    separated_conservative = (lam1 - delta1) > (lam2 + delta2)

    return {
        "N": N,
        "lambda1": lam1,
        "lambda2": lam2,
        "delta1": delta1,
        "delta2": delta2,
        "simple_test": separated_simple,
        "conservative_test": separated_conservative}

def effective_N(pc):
    r1 = np.corrcoef(pc[:-1], pc[1:])[0,1]
    N = len(pc)
    return N * (1 - r1) / (1 + r1)

res_E1 = north_test(VAR_E1, PC_E1[0,:], mode1=0, mode2=1)
res_E2 = north_test(VAR_E2, PC_E2[0,:], mode1=0, mode2=1)

print("=== AMOC ON: EOF1 vs EOF2 ===")
for k, v in res_E1.items():
    print(f"{k}: {v}")

print("\n=== AMOC OFF: EOF1 vs EOF2 ===")
for k, v in res_E2.items():
    print(f"{k}: {v}")

# %%

plt.figure()
plt.plot(PC_E1[0,:], label='First PC')
plt.plot(PC_E1[1,:], label='Second PC')
plt.xlabel('Time [model year]')
plt.ylabel('PC')
plt.title('AMOC on ')
plt.legend()

plt.figure()
plt.plot(PC_E2[0,:], label='First PC')
plt.plot(PC_E2[1,:], label='Second PC')
plt.xlabel('Time [model year]')
plt.ylabel('PC')
plt.title('AMOC off')
plt.legend()

# %% SST variability

fh = netcdf.Dataset(directory + 'SST_SSS_year_1_100_SO.nc','r')

time_E1          = fh.variables['time'][:]         #time [model years]
lon             = fh.variables['lon'][:]          #longitude [degrees_east]
lat             = fh.variables['lat'][:]          #latitude [degrees_north]
SST_E1           = fh.variables['SST'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'SST_SSS_year_500_600_SO.nc','r')

time_E2          = fh.variables['time'][:]         #time [model years]
SST_E2           = fh.variables['SST'][:]    #Transport [Sv]

fh.close()

#%%

plt.figure()
plt.contourf(lon, lat, np.std(SST_E1[:,:,:], axis=0))

plt.figure()
plt.contourf(np.std(SST_E2[:,:,:], axis=0))


# %%

import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Subsample the data
lon_sub = lon[::5]  # Take every 5th longitude point
lat_sub = lat[::5]  # Take every 5th latitude point
SST_SO_1_sub = SST_E1[:, ::5, ::5]  # Subsample the data
SST_SO_3_sub = SST_E2[:, ::5, ::5]  # Subsample the data

# Create a 1x2 subplot with South Polar Stereographic projection
fig, axs = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.SouthPolarStereo()})

# Add coastlines and set extent
for ax in axs:
    ax.coastlines(resolution='110m')  # Use lower resolution
    ax.set_extent([-180, 180, -90, -30], crs=ccrs.PlateCarree())  # Focus on Antarctica
    
gl = axs[0].gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

gl = axs[1].gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

axs[0].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none'
)

axs[1].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none'
)

c1 = axs[0].pcolormesh(
    lon_sub, lat_sub, np.var(SST_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-1, vmax=1
)
fig.colorbar(c1, ax=axs[0], orientation='vertical', shrink=0.8)
axs[0].set_title('a) SST variance - AMOC on', fontsize=16)

c2 = axs[1].pcolormesh(
    lon_sub, lat_sub, np.var(SST_SO_3_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-1, vmax=1
)
fig.colorbar(c2, ax=axs[1], orientation='vertical', shrink=0.8)
axs[1].set_title('b) SST variance - AMOC off', fontsize=16)

#Define region for som regions
lon_min, lon_max = -50, 0
lat_min, lat_max = -50, -35

lon_min_NZ, lon_max_NZ = 150, 190
lat_min_NZ, lat_max_NZ = -70, -60

lon_min_AU, lon_max_AU = 80, 150
lat_min_AU, lat_max_AU = -70, -40

lon_min_PA, lon_max_PA = -110, -60
lat_min_PA, lat_max_PA = -70, -40

lon_min_p, lon_max_p = 170, 250
lat_min_p, lat_max_p = -60, -45

#ACC transect
AT_lon_min, AT_lon_max = -60, 25
AT_lat_min, AT_lat_max = -90, -30

PA_lon_min = 150

x_1_AT  = np.arange(AT_lon_min, AT_lon_max + 0.1, 0.1)
y_1_AT	= np.zeros(len(x_1_AT)) + AT_lat_min
y_2_AT	= np.zeros(len(x_1_AT)) + AT_lat_max

y_3_AT = np.arange(AT_lat_min, AT_lat_max + 0.1, 0.1)
x_2_AT = np.zeros(len(y_3_AT)) + AT_lon_min
x_3_AT = np.zeros(len(y_3_AT)) + AT_lon_max
x_1_PA = np.zeros(len(y_3_AT)) + PA_lon_min

#WGKP
lon_min_wgkp, lon_max_wgkp = -35, 80
lat_min_wgkp, lat_max_wgkp = -80, -50

x_1	= np.arange(lon_min, lon_max + 0.1, 0.1)
y_1	= np.zeros(len(x_1)) + lat_min
y_2	= np.zeros(len(x_1)) + lat_max

y_3 = np.arange(lat_min, lat_max + 0.1, 0.1)
x_2 = np.zeros(len(y_3)) + lon_min
x_3 = np.zeros(len(y_3)) + lon_max

axs[0].plot(x_1, y_1, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_1, y_2, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_2, y_3, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_3, y_3, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[1].plot(x_1, y_1, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_1, y_2, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_2, y_3, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_3, y_3, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_p	= np.arange(lon_min_p, lon_max_p + 0.1, 0.1)
y_1_p	= np.zeros(len(x_1_p)) + lat_min_p
y_2_p	= np.zeros(len(x_1_p)) + lat_max_p

y_3_p = np.arange(lat_min_p, lat_max_p + 0.1, 0.1)
x_2_p = np.zeros(len(y_3_p)) + lon_min_p
x_3_p = np.zeros(len(y_3_p)) + lon_max_p

axs[0].plot(x_1_p, y_1_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_1_p, y_2_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_2_p, y_3_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_3_p, y_3_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[1].plot(x_1_p, y_1_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_1_p, y_2_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_2_p, y_3_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_3_p, y_3_p, '--k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'SST_variance_diff_contours_regions.pdf')
plt.show()

# %%
