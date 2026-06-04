#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  2 16:12:53 2025

@author: 6008399

MLD, SST and SSS in Drake Passage region 

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
from scipy.optimize import fsolve
from numpy.polynomial.polynomial import Polynomial
from scipy.io import loadmat
from statsmodels.graphics.gofplots import qqplot

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/HR_POP/netcdf/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'


#%%

fh = netcdf.Dataset(directory + 'MLD_year_1-600_Australia.nc','r')

time            = fh.variables['time'][:]       #time [model years]
MLD_AU             = fh.variables['MLD'][:]        #Transport [Sv]
MLD_max_AU         = fh.variables['MLD_max'][:]
lon_AU             = fh.variables['lon'][:]        #Longitude
lat_AU             = fh.variables['lat'][:]        #Latitude

fh.close()

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_wgkp= fh.variables['MLD_max'][:] 

fh.close()

#%%

fig, axs = plt.subplots(1, 3, figsize=(14, 4))  

contour1 = axs[0].contourf(lon, lat, np.mean(MLD_AU[63:114, :, :], axis=0), levels=np.linspace(0,600,11), cmap='viridis')
fig.colorbar(contour1, ax=axs[0], orientation='vertical', label='Mean MLD')
axs[0].set_title('SOM cycle 1')
axs[0].set_xlabel('Longitude [$^\circ$E]')
axs[0].set_ylabel('Latitude [$^\circ$N]')

contour2 = axs[1].contourf(lon, lat, np.mean(MLD_AU[324:378, :, :], axis=0), levels=np.linspace(0,600,11), cmap='viridis')
fig.colorbar(contour2, ax=axs[1], orientation='vertical', label='Mean MLD')
axs[1].set_title('SOM cycle 2')
axs[1].set_xlabel('Longitude [$^\circ$E]')
#axs[1].set_ylabel('Latitude [$^\circ$N]')

contour3 = axs[2].contourf(lon, lat, np.mean(MLD_AU[500:600, :, :], axis=0), levels=np.linspace(0,600,11), cmap='viridis')
fig.colorbar(contour3, ax=axs[2], orientation='vertical', label='Mean MLD')
axs[2].set_title('SOM cycle 3')
axs[2].set_xlabel('Longitude [$^\circ$E]')
#axs[2].set_ylabel('Latitude [$^\circ$N]')

fig.suptitle('Mean MLD', fontsize=15)

plt.tight_layout()
plt.show()

#%%

fig, axs = plt.subplots(1, 3, figsize=(14, 4))  

contour1 = axs[0].contourf(lon, lat, np.mean(MLD[324:378, :, :], axis=0) - np.mean(MLD[63:114, :, :], axis=0), levels=np.linspace(-350,350,19), cmap='seismic')
fig.colorbar(contour1, ax=axs[0], orientation='vertical', label='Mean MLD difference')
axs[0].set_title('SOM cycle 2 - 1')
axs[0].set_xlabel('Longitude [$^\circ$E]')
axs[0].set_ylabel('Latitude [$^\circ$N]')

contour2 = axs[1].contourf(lon, lat, np.mean(MLD[500:600, :, :], axis=0) - np.mean(MLD[324:378, :, :], axis=0), levels=np.linspace(-350,350,19), cmap='seismic')
fig.colorbar(contour2, ax=axs[1], orientation='vertical', label='Mean MLD difference')
axs[1].set_title('SOM cycle 3 - 2')
axs[1].set_xlabel('Longitude [$^\circ$E]')
#axs[1].set_ylabel('Latitude [$^\circ$N]')

contour3 = axs[2].contourf(lon, lat, np.mean(MLD[500:600, :, :], axis=0) - np.mean(MLD[63:114, :, :], axis=0), levels=np.linspace(-350,350,19), cmap='seismic')
fig.colorbar(contour3, ax=axs[2], orientation='vertical', label='Mean MLD difference')
axs[2].set_title('SOM cycle 3 - 1')
axs[2].set_xlabel('Longitude [$^\circ$E]')
#axs[2].set_ylabel('Latitude [$^\circ$N]')

fig.suptitle('MLD difference', fontsize=15)


plt.tight_layout()
plt.savefig(directory_figures + 'mean_mld_difference_AU_HRPOP.pdf')
plt.show()

#%%

plt.figure()
plt.plot(time, MLD_max, label='Max MLD (80 - 150E, 70 - 50S)')
plt.legend()

#%%

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/SOM-P.nc','r')

time        = fh.variables['time'][:] #model years
som_p         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh1 = netcdf.Dataset(directory + 'AMOC_transport_depth_0-1000m.nc','r')

amoc_strength   = fh1.variables['Transport'][:]     #Volume transport [Sv]

fh1.close()

fh = netcdf.Dataset(directory + 'Drake_Passage_transport.nc','r')

time                = fh.variables['time'][:]         #time [model years]
acc_transport           = fh.variables['Transport'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_depth_0-1000m.nc','r')

time                = fh.variables['time'][:]         #time [model years]
som_star           = fh.variables['SOM'][:]    #Transport [Sv]

fh.close()

#%%

#Central moving average
def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

fig, ax1 = plt.subplots(figsize=(8, 5))

# Plot the first variable (MLD_max) on the primary y-axis
ax1.plot(Moving_average(MLD_max, 20), color='blue', label='Max MLD (80 - 150E, 70 - 50S)')
ax1.set_xlabel('Time [years]', fontsize=12)
ax1.set_ylabel('Max MLD [units]', color='blue', fontsize=12)
ax1.tick_params(axis='y', labelcolor='blue')
ax1.legend(loc='upper left')

# Create the second y-axis
ax2 = ax1.twinx()  # Share the same x-axis
ax2.plot(time, acc_transport, color='red', label='ACC transport')
ax2.set_ylabel('ACC Transport [units]', color='red', fontsize=12)
ax2.tick_params(axis='y', labelcolor='red')
ax2.legend(loc='upper right')

# Add a title and grid
plt.title('Max MLD (moving average = 20) and ACC Transport', fontsize=14)
plt.grid()

# Show the plot
plt.tight_layout()
plt.savefig(directory_figures + 'max_mld_AU_acc_transport_HRPOP.pdf')
plt.show()

#%%

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data


plt.figure(figsize=(8,4))
plt.plot(time, norm(som, som), color='blue', label='SOM index')
plt.plot(time[window//2:-window//2 + 1], norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), color='black', label='max MLD')
plt.plot(time, norm(acc_transport, acc_transport), color='green', label='Drake Passage')
plt.grid()
plt.legend()
plt.xlabel('Time [model years]')
plt.ylabel('Normalized scale')
plt.title('Max MLD (80 - 150E, 70 - 50S)')
plt.savefig(directory_figures + 'ACC_SOM_max_MLD_Australie.pdf')


#%%

fh = netcdf.Dataset(directory + 'MLD_year_1-600_SO.nc','r')

time_1            = fh.variables['time'][0:100]       #time [model years]
MLD_SO_1             = fh.variables['MLD'][0:100, :,:]        #Transport [Sv]
time_3           = fh.variables['time'][500:600]       #time [model years]
MLD_SO_3             = fh.variables['MLD'][500:600, :,:]        #Transport [Sv]
MLD_max         = fh.variables['MLD_max'][:]
lon             = fh.variables['lon'][:]        #Longitude
lat             = fh.variables['lat'][:]        #Latitude
MLD_SO_NZ       = fh.variables['MLD'][:, 200:400,2600:3000]  
lon_NZ          = np.linspace(150, 190, len(lon[2600:3000]))
lat_NZ          = lat[200:400]
MLD_SO_PA       = fh.variables['MLD'][:, 200:600,0:500]  
lon_PA          = lon[0:500]
lat_PA          = lat[200:600]

fh.close()

#%%

plt.figure(figsize=(8,4))
plt.plot(time, norm(som, som), color='blue', label='SOM index')
plt.plot(time[window//2:-window//2 + 1], norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), color='black', label='max MLD')
plt.plot(time, norm(acc_transport, acc_transport), color='green', label='Drake Passage')
plt.grid()
plt.legend()
plt.xlabel('Time [model years]')
plt.ylabel('Normalized scale')
#plt.title('Max MLD (80 - 150E, 70 - 50S)')

#%%

fig, axs = plt.subplots(3, 1, figsize=(6, 8))  

contour1 = axs[0].contourf(lon, lat, np.mean(MLD_SO_1, axis=0), levels=np.linspace(0,500,19))
fig.colorbar(contour1, ax=axs[0], orientation='vertical', label='Mean MLD')
axs[0].set_title('a) Model year 1-100')
axs[0].set_xlabel('Longitude [$^\circ$E]')
axs[0].set_ylabel('Latitude [$^\circ$N]')

contour2 = axs[1].contourf(np.mean(MLD_SO_3, axis=0), levels=np.linspace(0,500,19))
fig.colorbar(contour2, ax=axs[1], orientation='vertical', label='Mean MLD')
axs[1].set_title('b) Model year 500-600')
axs[1].set_xlabel('Longitude [$^\circ$E]')
#axs[1].set_ylabel('Latitude [$^\circ$N]')

contour3 = axs[2].contourf(np.mean(MLD_SO_3, axis=0) - np.mean(MLD_SO_1, axis=0), levels=np.linspace(-350,350,19), cmap='seismic')
fig.colorbar(contour3, ax=axs[2], orientation='vertical', label='Mean MLD difference')
axs[2].set_title('c) Difference')
axs[2].set_xlabel('Longitude [$^\circ$E]')
#axs[2].set_ylabel('Latitude [$^\circ$N]')

fig.suptitle('MLD SO difference', fontsize=15)


plt.tight_layout()
#plt.savefig(directory_figures + 'mean_mld_difference_SO_HRPOP.pdf')
plt.show()

#%%

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# Subsample the data
lon_sub = lon[::5]  # Take every 5th longitude point
lat_sub = lat[::5]  # Take every 5th latitude point
MLD_SO_3_sub = MLD_SO_3[:, ::5, ::5]  # Subsample the data
MLD_SO_1_sub = MLD_SO_1[:, ::5, ::5]  # Subsample the data

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

# First subplot: Maximum MLD difference
c1 = axs[0].pcolormesh(
    lon_sub, lat_sub, np.max(MLD_SO_3_sub, axis=0) - np.max(MLD_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-300, vmax=300
)
fig.colorbar(c1, ax=axs[0], orientation='vertical', shrink=0.8)
axs[0].set_title('a) Maximum MLD difference', fontsize=16)

# Second subplot: MLD variance difference
c2 = axs[1].pcolormesh(
    lon_sub, lat_sub, np.var(MLD_SO_3_sub, axis=0) - np.var(MLD_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-5000, vmax=5000
)
fig.colorbar(c2, ax=axs[1], orientation='vertical', shrink=0.8)
axs[1].set_title('b) MLD variance difference', fontsize=16)

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

#axs[0].plot(x_1_AT, y_1_AT, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_1_PA, y_3_AT, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_2_AT, y_3_AT, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0].plot(x_3_AT, y_3_AT, '-k', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_AU	= np.arange(lon_min_AU, lon_max_AU + 0.1, 0.1)
y_1_AU	= np.zeros(len(x_1_AU)) + lat_min_AU
y_2_AU	= np.zeros(len(x_1_AU)) + lat_max_AU

y_3_AU = np.arange(lat_min_AU, lat_max_AU + 0.1, 0.1)
x_2_AU = np.zeros(len(y_3_AU)) + lon_min_AU
x_3_AU = np.zeros(len(y_3_AU)) + lon_max_AU

x_1_NZ	= np.arange(lon_min_NZ, lon_max_NZ + 0.1, 0.1)
y_1_NZ	= np.zeros(len(x_1_NZ)) + lat_min_NZ
y_2_NZ	= np.zeros(len(x_1_NZ)) + lat_max_NZ

y_3_NZ = np.arange(lat_min_NZ, lat_max_NZ + 0.1, 0.1)
x_2_NZ = np.zeros(len(y_3_NZ)) + lon_min_NZ
x_3_NZ = np.zeros(len(y_3_NZ)) + lon_max_NZ

x_1_PA	= np.arange(lon_min_PA, lon_max_PA + 0.1, 0.1)
y_1_PA	= np.zeros(len(x_1_PA)) + lat_min_PA
y_2_PA	= np.zeros(len(x_1_PA)) + lat_max_PA

y_3_PA = np.arange(lat_min_PA, lat_max_PA + 0.1, 0.1)
x_2_PA = np.zeros(len(y_3_PA)) + lon_min_PA
x_3_PA = np.zeros(len(y_3_PA)) + lon_max_PA

axs[1].plot(x_1_NZ, y_1_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_1_NZ, y_2_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_2_NZ, y_3_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_3_NZ, y_3_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[1].plot(x_1_AU, y_1_AU, '-', color='purple', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_1_AU, y_2_AU, '-', color='purple', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_2_AU, y_3_AU, '-', color='purple', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_3_AU, y_3_AU, '-', color='purple', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[1].plot(x_1_PA, y_1_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_1_PA, y_2_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_2_PA, y_3_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_3_PA, y_3_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

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

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp

axs[1].plot(x_1_wgkp, y_1_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_1_wgkp, y_2_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_2_wgkp, y_3_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[1].plot(x_3_wgkp, y_3_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'MLD_max_variance_diff_contours_regions.pdf')
plt.show()

#%%

# Subsample the data
lon_sub = lon[::5]  # Take every 5th longitude point
lat_sub = lat[::5]  # Take every 5th latitude point
MLD_SO_3_sub = MLD_SO_3[:, ::5, ::5]  # Subsample the data
MLD_SO_1_sub = MLD_SO_1[:, ::5, ::5]  # Subsample the data

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

# First subplot: Maximum MLD difference
c1 = axs[0].pcolormesh(
    lon_sub, lat_sub, np.max(MLD_SO_3_sub, axis=0) - np.max(MLD_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-300, vmax=300
)
fig.colorbar(c1, ax=axs[0], orientation='vertical', shrink=0.8)
axs[0].set_title('a) Maximum MLD difference', fontsize=16)

# Second subplot: MLD variance difference
c2 = axs[1].pcolormesh(
    lon_sub, lat_sub, np.var(MLD_SO_3_sub, axis=0) - np.var(MLD_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-5000, vmax=5000
)
fig.colorbar(c2, ax=axs[1], orientation='vertical', shrink=0.8)
axs[1].set_title('b) MLD variance difference', fontsize=16)

# Adjust the layout
plt.tight_layout()
plt.savefig(directory_figures +'MLD_max_variance_diff.pdf')
plt.show()
#%%

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

# Example data
lon = np.linspace(0, 360, 100)
lat = np.linspace(-90, -35, 50)
MLD_SO_3 = np.random.rand(12, 50, 100) * 100
MLD_SO_1 = np.random.rand(12, 50, 100) * 100

# Create a 1x2 subplot with South Polar Stereographic projection
fig, axs = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': ccrs.SouthPolarStereo()})

# Add coastlines and set extent
for ax in axs:
    ax.coastlines(resolution='110m')  # Use lower resolution
    ax.set_extent([-180, 180, -90, -60], crs=ccrs.PlateCarree())  # Focus on Antarctica

# Add gridlines
for ax in axs:
    gl = ax.gridlines(draw_labels=True, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = False
    gl.xlabel_style = {'size': 10, 'color': 'black'}
    gl.ylabel_style = {'size': 10, 'color': 'black'}

# First subplot: Maximum MLD difference
c1 = axs[0].contourf(
    lon, lat, np.max(MLD_SO_3, axis=0) - np.max(MLD_SO_1, axis=0),
    transform=ccrs.PlateCarree(), levels=np.linspace(-350, 350, 19), cmap='RdBu_r', extend='both'
)
fig.colorbar(c1, ax=axs[0], orientation='vertical', shrink=1)
axs[0].set_title('a) Maximum MLD difference', fontsize=13)

# Second subplot: MLD variance difference
c2 = axs[1].contourf(
    lon, lat, np.var(MLD_SO_3, axis=0) - np.var(MLD_SO_1, axis=0),
    transform=ccrs.PlateCarree(), levels=np.linspace(-5000, 5000, 19), cmap='RdBu_r', extend='both'
)
fig.colorbar(c2, ax=axs[1], orientation='vertical', shrink=1)
axs[1].set_title('b) MLD variance difference', fontsize=13)

# Adjust layout
plt.tight_layout()
plt.show()





#%%

fig, axs = plt.subplots(3, 1, figsize=(6, 8))  

contour1 = axs[0].contourf(np.var(MLD_SO_1, axis=0), levels=np.linspace(0,8000,19), extend = 'both')
fig.colorbar(contour1, ax=axs[0], orientation='vertical', label='MLD variance [m$^2$]')
axs[0].set_title('a) Model year 1-100')
axs[0].set_xlabel('Longitude [$^\circ$E]')
axs[0].set_ylabel('Latitude [$^\circ$N]')

contour2 = axs[1].contourf(np.var(MLD_SO_3, axis=0), levels=np.linspace(0,8000,19), extend='both')
fig.colorbar(contour2, ax=axs[1], orientation='vertical', label='MLD variance [m$^2$]')
axs[1].set_title('b) Model year 500-600')
axs[1].set_xlabel('Longitude [$^\circ$E]')
#axs[1].set_ylabel('Latitude [$^\circ$N]')

contour3 = axs[2].contourf(np.var(MLD_SO_3, axis=0) - np.var(MLD_SO_1, axis=0), levels=np.linspace(-5000,5000,19), extend='both', cmap='seismic')
fig.colorbar(contour3, ax=axs[2], orientation='vertical', label='MLD variance difference [m$^2$]')
axs[2].set_title('c) Difference')
axs[2].set_xlabel('Longitude [$^\circ$E]')
#axs[2].set_ylabel('Latitude [$^\circ$N]')

fig.suptitle('MLD variance SO difference', fontsize=15)


plt.tight_layout()
plt.savefig(directory_figures + 'variance_mld_difference_SO_HRPOP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(3, 1, figsize=(6, 8))  

contour1 = axs[0].contourf(lon_NZ, lat_NZ, np.var(MLD_SO_NZ[500:600], axis=0) - np.var(MLD_SO_NZ[0:100], axis=0), levels=np.linspace(-5000,5000,19), extend='both', cmap='seismic')
fig.colorbar(contour1, ax=axs[0], orientation='vertical', label='MLD variance difference [m$^2$]')
axs[0].set_title('a) NZ')
axs[0].set_xlabel('Longitude [$^\circ$E]')
axs[0].set_ylabel('Latitude [$^\circ$N]')

contour2 = axs[1].contourf(lon_AU, lat_AU, np.var(MLD_AU[500:600], axis=0) - np.var(MLD_AU[0:100], axis=0), levels=np.linspace(-5000,5000,19), extend='both', cmap='seismic')
fig.colorbar(contour2, ax=axs[1], orientation='vertical', label='MLD variance difference [m$^2$]')
axs[1].set_title('b) AU')
axs[1].set_xlabel('Longitude [$^\circ$E]')
#axs[1].set_ylabel('Latitude [$^\circ$N]')

contour3 = axs[2].contourf(lon_PA, lat_PA, np.var(MLD_SO_PA[500:600], axis=0) - np.var(MLD_SO_PA[0:100], axis=0), levels=np.linspace(-5000,5000,19), extend='both', cmap='seismic')
fig.colorbar(contour3, ax=axs[2], orientation='vertical', label='MLD variance difference [m$^2$]')
axs[2].set_title('c) Pacific')
axs[2].set_xlabel('Longitude [$^\circ$E]')
#axs[2].set_ylabel('Latitude [$^\circ$N]')

fig.suptitle('MLD variance difference', fontsize=15)


plt.tight_layout()
plt.savefig(directory_figures + 'variance_mld_difference_AU_NZ_PA_HRPOP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(3, 1, figsize=(6, 8))  

contour1 = axs[0].contourf(np.max(MLD_SO_1, axis=0), levels=np.linspace(0,500,19))
fig.colorbar(contour1, ax=axs[0], orientation='vertical', label='Max MLD')
axs[0].set_title('a) Model year 1-100')
axs[0].set_xlabel('Longitude [$^\circ$E]')
axs[0].set_ylabel('Latitude [$^\circ$N]')

contour2 = axs[1].contourf(np.max(MLD_SO_3, axis=0), levels=np.linspace(0,500,19))
fig.colorbar(contour2, ax=axs[1], orientation='vertical', label='Max MLD')
axs[1].set_title('b) Model year 500-600')
axs[1].set_xlabel('Longitude [$^\circ$E]')
#axs[1].set_ylabel('Latitude [$^\circ$N]')

contour3 = axs[2].contourf(np.max(MLD_SO_3, axis=0) - np.max(MLD_SO_1, axis=0), levels=np.linspace(-350,350,19), cmap='seismic')
fig.colorbar(contour3, ax=axs[2], orientation='vertical', label='Max MLD difference')
axs[2].set_title('c) Difference')
axs[2].set_xlabel('Longitude [$^\circ$E]')
#axs[2].set_ylabel('Latitude [$^\circ$N]')

fig.suptitle('MLD SO difference', fontsize=15)


plt.tight_layout()
plt.savefig(directory_figures + 'max_mld_difference_SO_HRPOP.pdf')
plt.show()

#%%

plt.figure()

plt.plot(MLD_max_AU, label='AU')
#plt.plot(MLD_max_wgkp, label='WGKP')
#plt.plot(MLD_max, label='SO')
#plt.plot(np.max(MLD_SO_NZ, axis=(1,2)), label='NZ')
plt.plot(np.max(MLD_SO_PA, axis=(1,2)), label='PA')
plt.ylim(2000,0)
plt.title('Max MLD')
plt.legend()

#%%

plt.figure()
plt.plot(np.mean(MLD, axis=(1,2)), label='AU')
plt.plot(np.mean(MLD_SO_NZ, axis=(1,2)), label='NZ')
plt.plot(np.mean(MLD_SO_PA, axis=(1,2)), label='PA')
plt.ylim(300,0)
plt.title('Mean MLD')
plt.legend()

#%%

plt.figure()
plt.plot(np.mean(MLD_SO_NZ, axis=(1,2)), label='mean')
plt.plot(Moving_average(np.max(MLD_SO_NZ, axis=(1,2)), 20), label='max')
plt.plot(np.min(MLD_SO_NZ, axis=(1,2)), label='min')
plt.ylim(1100,0)
plt.title('MLD NZ')
plt.legend()

plt.figure()
plt.plot(np.mean(MLD_AU, axis=(1,2)), label='mean')
plt.plot(Moving_average(np.max(MLD_AU, axis=(1,2)), 20), label='max')
plt.plot(np.min(MLD_AU, axis=(1,2)), label='min')
plt.ylim(600,0)
plt.title('MLD AU')
plt.legend()

plt.figure()
plt.plot(np.mean(MLD_SO_PA, axis=(1,2)), label='mean')
plt.plot(Moving_average(np.max(MLD_SO_PA, axis=(1,2)), 20), label='max')
plt.plot(np.min(MLD_SO_PA, axis=(1,2)), label='min')
plt.ylim(600,0)
plt.title('MLD PA')
plt.legend()

fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# ---- 1) New Zealand sector ----
axs[0].plot(np.mean(MLD_SO_NZ, axis=(1,2)), label='mean')
axs[0].plot(Moving_average(np.max(MLD_SO_NZ, axis=(1,2)), 20), label='max')
#axs[0].plot(np.min(MLD_SO_NZ, axis=(1,2)), label='min')
axs[0].set_ylim(1100, 0)
axs[0].set_title('MLD NZ')
axs[0].legend()

# ---- 2) Australian sector ----
axs[1].plot(np.mean(MLD_AU, axis=(1,2)), label='mean')
axs[1].plot(Moving_average(np.max(MLD_AU, axis=(1,2)), 20), label='max')
#axs[1].plot(np.min(MLD_AU, axis=(1,2)), label='min')
axs[1].set_ylim(600, 0)
axs[1].set_title('MLD AU')
axs[1].legend()

# ---- 3) Pacific sector ----
axs[2].plot(np.mean(MLD_SO_PA, axis=(1,2)), label='mean')
axs[2].plot(Moving_average(np.max(MLD_SO_PA, axis=(1,2)), 20), label='max')
#axs[2].plot(np.min(MLD_SO_PA, axis=(1,2)), label='min')
axs[2].set_ylim(600, 0)
axs[2].set_title('MLD PA')
axs[2].legend()

plt.tight_layout()
plt.show()


#%%

# Find the maximum value in the array
max_value = np.max(MLD_SO_3[-1])

# Find the x and y coordinates of the maximum value
x, y = np.where(MLD_SO_3[-1] == max_value)

# Print the results
print(f"Maximum value: {max_value}")
print(f"Coordinates of maximum value: x={x[0]}, y={y[0]}")

#%%


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.animation import FFMpegWriter

# Example data
time = time_3
#lon[lon<0]	= lon[lon<0] + 360.0
lon_grid, lat_grid = np.meshgrid(lon, lat)
data = MLD_SO_3
line_data_1 = np.max(MLD_SO_PA[500::], axis=(1,2))
line_data_2 = MLD_max_AU[500::]
line_data_3 = MLD_max_NZ[500::]

# Create the figure and subplots
fig = plt.figure(figsize=(8, 6))

# Subplot 1: Contour plot
ax1 = fig.add_subplot(2, 1, 1)#, projection=ccrs.PlateCarree())
#ax1 = fig.add_subplot(1, 1, 1)
contour = ax1.contourf(data[0], levels=np.linspace(0, 600, 31), extend='both')#, transform=ccrs.PlateCarree())
#ax1.coastlines()
#ax1.add_feature(cfeature.LAND)
#ax1.set_title('Contour Plot')
cbar = plt.colorbar(contour, ax=ax1, orientation='vertical', pad=0.05)
cbar.set_label('MLD')

ax2 = fig.add_subplot(2, 1, 2)
line_1, = ax2.plot(time, line_data_1, color='blue', label='PA') 
line_2, = ax2.plot(time, line_data_2, color='orange', label='AU') 
line_3, = ax2.plot(time, line_data_3, color='green', label='NZ') 
ax2.set_xlim(time[0], time[-1])
#ax2.set_ylim(1.22e18, 1.35e18)
ax2.set_title('MLD maximum') 
ax2.set_xlabel('Time [model year]')
ax2.set_ylabel('MLD [m]')
ax2.legend()

#%%

# Update function for animation
def update(frame):
    # Update contour plot
    ax1.clear()
    contour = ax1.contourf(data[frame], levels=np.linspace(0, 600, 31),extend = 'both')
    #ax1.coastlines()
    #ax1.add_feature(cfeature.LAND)
    ax1.set_title(f'MLD SOM cycle 3 (Time: {frame} model year)')

    # Update line plot
    line_1.set_data(time[:frame], line_data_1[:frame])
    line_2.set_data(time[:frame], line_data_2[:frame])
    line_3.set_data(time[:frame], line_data_3[:frame])
    return contour, line

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=len(time), interval=100)
ani.save(directory_figures + 'MLD_SOM_3_SO.gif', writer='ffmpeg', fps=5)

# Show the animation
plt.show() 


