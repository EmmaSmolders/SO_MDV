#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 18 14:34:50 2025

@author: 6008399

Comparison of densities in the different convection zones


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
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/HR_POP/netcdf/'
directory_zenodo = r'/Users/6008399/Documents/PhD/HR_POP/Zenodo/Data_Final/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%% Functions

def compute_N2_from_profile(density, depth, rho0=1027, g=9.81):
    """
    Compute vertical density gradient (drho/dz) and buoyancy frequency N²
    from a 1D density profile and depth array.

    Compute N² exactly following the user's original code,
    assuming depth increases downward.
    """

    PA = len(depth)
    n0 = np.zeros(PA)

    for i in range(PA):
        if i == 0:
            # surface: (rho0 - rho1)/(z1 - z0)
            n0[i] = (density[i] - density[i+1]) / (depth[i+1] - depth[i])

        elif i == PA - 1:
            # bottom: (rho_{n-2} - rho_{n-1})/(z_{n-1}-z_{n-2})
            n0[i] = (density[i-1] - density[i]) / (depth[i] - depth[i-1])

        else:
            # centered: (rho[i-1] - rho[i+1]) / (z[i+1] - z[i-1])
            n0[i] = (density[i-1] - density[i+1]) / (depth[i+1] - depth[i-1])

    # Your formula: N² = -(g / rho0) * n0
    N2 = - (g / rho0) * n0
    
    if N2.any() < 0:
        print('N2 is unstable!')

    return n0, N2

def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#%% Read in data 

fh = netcdf.Dataset(directory_zenodo + 'PD_gsw_year_1-600_area_averaged_NZ_70_62S_160E_190E.nc','r')

depth_PD       = fh.variables['depth'][:]  #depth
time_all = fh.variables['time'][:] #Time SOM cycle 1
PD_NZ_all = fh.variables['PD'][:]
PD_NZ_1 = fh.variables['PD'][1:100]
PD_NZ_2 = fh.variables['PD'][500:600]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/PD_gsw_year_1-600_area_averaged_WGKP.nc','r')

depth_PD       = fh.variables['depth'][:]  #depth
time_all = fh.variables['time'][:] #Time SOM cycle 1
PD_WGKP_all = fh.variables['PD'][:]
PD_WGKP_1 = fh.variables['PD'][1:100]
PD_WGKP_2 = fh.variables['PD'][500:600]

fh.close()

fh = netcdf.Dataset(directory_zenodo + 'PD_gsw_year_1-600_area_averaged_AU_70_50S_80_150E.nc','r')
depth_PD        = fh.variables['depth'][:]  #depth
time_all        = fh.variables['time'][:] #Time SOM cycle 1
PD_AU_all       = fh.variables['PD'][:]
PD_AU_1       = fh.variables['PD'][1:100]
PD_AU_2       = fh.variables['PD'][500:600]

fh.close()

fh = netcdf.Dataset(directory_zenodo + 'PD_gsw_year_1-600_area_averaged_PA_70_50S_110W_60W.nc','r')

depth_PD        = fh.variables['depth'][:]  #depth
time_all        = fh.variables['time'][:] #Time SOM cycle 1
PD_PA_all       = fh.variables['PD'][:]
PD_PA_1       = fh.variables['PD'][1:100]
PD_PA_2       = fh.variables['PD'][500:600]

fh.close()

fh = netcdf.Dataset(directory + 'MLD_year_1-600_SO.nc','r')

time_1            = fh.variables['time'][0:100]       #time [model years]
MLD_SO_1             = fh.variables['MLD'][0:100, :,:]        #Transport [Sv]
time_3           = fh.variables['time'][500:600]       #time [model years]
MLD_SO_3             = fh.variables['MLD'][500:600, :,:]        #Transport [Sv]
MLD_max         = fh.variables['MLD_max'][:]
lon             = fh.variables['lon'][:]        #Longitude
lat             = fh.variables['lat'][:]        #Latitude

fh.close()

#%% Determine N^2

n0_som_1_WGKP, N2_som_1_WGKP = compute_N2_from_profile(np.mean(PD_WGKP_1, axis = 0), depth)
n0_som_2_WGKP, N2_som_2_WGKP = compute_N2_from_profile(np.mean(PD_WGKP_2, axis=0), depth)

n0_som_1_NZ, N2_som_1_NZ = compute_N2_from_profile(np.mean(PD_NZ_1, axis = 0), depth)
n0_som_2_NZ, N2_som_2_NZ = compute_N2_from_profile(np.mean(PD_NZ_2, axis=0), depth)

n0_som_1_AU, N2_som_1_AU = compute_N2_from_profile(np.mean(PD_AU_1, axis = 0), depth)
n0_som_2_AU, N2_som_2_AU = compute_N2_from_profile(np.mean(PD_AU_2, axis=0), depth)

n0_som_1_PA, N2_som_1_PA = compute_N2_from_profile(np.mean(PD_PA_1, axis = 0), depth)
n0_som_2_PA, N2_som_2_PA = compute_N2_from_profile(np.mean(PD_PA_2, axis=0), depth)


#%%

#Subsample the data otherwise takes a very long time to plot
lon_sub = lon[::5]  #Take every 5th longitude point etc.
lat_sub = lat[::5] 
MLD_SO_3_sub = MLD_SO_3[:, ::5, ::5]  
MLD_SO_1_sub = MLD_SO_1[:, ::5, ::5]  

fig = plt.figure(figsize=(10, 12))

axs = np.empty((3,2), dtype=object)

# Top row — South Polar maps
axs[0,0] = fig.add_subplot(3, 2, 1, projection=ccrs.SouthPolarStereo())
axs[0,1] = fig.add_subplot(3, 2, 2, projection=ccrs.SouthPolarStereo())

# Bottom row — normal (Cartesian) plots
axs[1,0] = fig.add_subplot(3, 2, 3)
axs[1,1] = fig.add_subplot(3, 2, 4)

axs[2,0] = fig.add_subplot(3, 2, 5)
axs[2,1] = fig.add_subplot(3, 2, 6)

# Add coastlines and gridlines to each subplot
axs[0,0].coastlines()
axs[0,0].set_extent([-180, 180, -90, -35], crs=ccrs.PlateCarree())  # Focus on Antarctica

axs[0,1].coastlines()
axs[0,1].set_extent([-180, 180, -90, -35], crs=ccrs.PlateCarree())  # Focus on Antarctica

    # Add gridlines
gl = axs[0,0].gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

gl = axs[0,1].gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

axs[0,0].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none'
)

axs[0,1].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none'
)

c1 = axs[0,0].pcolormesh(
    lon_sub, lat_sub, np.max(MLD_SO_3_sub, axis=0) - np.max(MLD_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-300, vmax=300)

cbar1 = fig.colorbar(c1, ax=axs[0, 0], orientation='vertical', extend='both', shrink=1)
cbar1.set_label('MLD maximum difference [m]', fontsize=12)
axs[0,0].set_title('a) MLD maximum difference', fontsize=14)

c2 = axs[0,1].pcolormesh(
    lon_sub, lat_sub, np.var(MLD_SO_3_sub, axis=0) - np.var(MLD_SO_1_sub, axis=0),transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-5000, vmax=5000)

cbar2 = fig.colorbar(c2, ax=axs[0,1], orientation='vertical', extend = 'both', shrink=1)
cbar2.set_label('MLD variance difference [m$^2$]', fontsize=12)
axs[0,1].set_title('b) MLD variance difference', fontsize=14)

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

#WGKP
lon_min_wgkp, lon_max_wgkp = -35, 80
lat_min_wgkp, lat_max_wgkp = -80, -50

x_1	= np.arange(lon_min, lon_max + 0.1, 0.1)
y_1	= np.zeros(len(x_1)) + lat_min
y_2	= np.zeros(len(x_1)) + lat_max

y_3 = np.arange(lat_min, lat_max + 0.1, 0.1)
x_2 = np.zeros(len(y_3)) + lon_min
x_3 = np.zeros(len(y_3)) + lon_max

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

axs[0,1].plot(x_1_NZ, y_1_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_NZ, y_2_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_NZ, y_3_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_NZ, y_3_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[0,1].plot(x_1_AU, y_1_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_AU, y_2_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_AU, y_3_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_AU, y_3_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[0,1].plot(x_1_PA, y_1_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_PA, y_2_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_PA, y_3_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_PA, y_3_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_p	= np.arange(lon_min_p, lon_max_p + 0.1, 0.1)
y_1_p	= np.zeros(len(x_1_p)) + lat_min_p
y_2_p	= np.zeros(len(x_1_p)) + lat_max_p

y_3_p = np.arange(lat_min_p, lat_max_p + 0.1, 0.1)
x_2_p = np.zeros(len(y_3_p)) + lon_min_p
x_3_p = np.zeros(len(y_3_p)) + lon_max_p

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp

axs[0,1].plot(x_1_wgkp, y_1_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_wgkp, y_2_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_wgkp, y_3_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_wgkp, y_3_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)


axs[1,0].plot(np.mean(PD_NZ_1, axis=0), depth, color='magenta', linewidth = 2.0, label='NZ')
axs[1,0].plot(np.mean(PD_WGKP_1, axis=0), depth, color='green', linewidth = 2.0, label='WGKP')
axs[1,0].plot(np.mean(PD_AU_1, axis=0), depth, color='darkred', linewidth = 2.0, label='AU')
axs[1,0].plot(np.mean(PD_PA_1, axis=0), depth, color='blue', linewidth = 2.0, label='PA')
axs[1,0].legend()
axs[1,0].set_ylim(2000, 1)  
axs[1,0].set_xlim(1026.4, 1027.8)
axs[1,0].set_title('c) Area-averaged PD (model year 1-100)', fontsize=14)
axs[1,0].set_xlabel('Potential density [kg/m$^3$]', fontsize=12)
axs[1,0].set_ylabel('Depth [m]', fontsize = 12)
axs[1,0].grid()

axs[1,1].plot(np.mean(PD_NZ_2, axis=0), depth, color='magenta', linewidth = 2.0, label='NZ')
axs[1,1].plot(np.mean(PD_WGKP_2, axis=0), depth, color='green', linewidth = 2.0, label='WGKP')
axs[1,1].plot(np.mean(PD_AU_2, axis=0), depth, color='darkred', linewidth = 2.0, label='AU')
axs[1,1].plot(np.mean(PD_PA_2, axis=0), depth, color='blue', linewidth = 2.0, label='PA')
axs[1,1].legend()
axs[1,1].set_ylim(2000, 1)  
axs[1,1].set_xlim(1026.4, 1027.8)
axs[1,1].set_title('d) Area-averaged PD (model year 500-600)', fontsize=14)
axs[1,1].set_xlabel('Potential density [kg/m$^3$]', fontsize=12)
axs[1,1].grid()

axs[2,0].plot(N2_som_1_NZ, depth, color='magenta', linewidth = 2.0, label='NZ')
axs[2,0].plot(N2_som_1_WGKP, depth, color='green', linewidth = 2.0, label='WGKP')
axs[2,0].plot(N2_som_1_AU, depth, color='darkred', linewidth = 2.0, label='AU')
axs[2,0].plot(N2_som_1_PA, depth, color='blue', linewidth = 2.0, label='PA')
axs[2,0].legend()
axs[2,0].set_ylim(2000, 1)  
#axs[2,0].set_xlim(1026.4, 1027.8)
axs[2,0].set_title('e) Area-averaged PD (model year 1-100)', fontsize=14)
axs[2,0].set_xlabel('$N^2$ [s$^{-1}$]', fontsize=12)
axs[2,0].set_ylabel('Depth [m]', fontsize = 12)
axs[2,0].grid()

axs[2,1].plot(N2_som_2_NZ, depth, color='magenta', linewidth = 2.0, label='NZ')
axs[2,1].plot(N2_som_2_WGKP, depth, color='green', linewidth = 2.0, label='WGKP')
axs[2,1].plot(N2_som_2_AU, depth, color='darkred', linewidth = 2.0, label='AU')
axs[2,1].plot(N2_som_2_PA, depth, color='blue', linewidth = 2.0, label='PA')
axs[2,1].legend()
axs[2,1].set_ylim(2000, 1)  
#axs[2,1].set_xlim(1026.4, 1027.8)
axs[2,1].set_title('f) Area-averaged $N^2$ (model year 500-600)', fontsize=14)
axs[2,1].set_xlabel('$N^2$ [s$^{-1}$]', fontsize=12)
axs[2,1].grid()

plt.tight_layout()
plt.savefig(directory_figures +'Figure_7_SOM_OS.pdf')
plt.show()


#%%

#Subsample the data otherwise takes a very long time to plot
lon_sub = lon[::5]  #Take every 5th longitude point etc.
lat_sub = lat[::5] 
MLD_SO_3_sub = MLD_SO_3[:, ::5, ::5]  
MLD_SO_1_sub = MLD_SO_1[:, ::5, ::5]  

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
gl = axs[0,0].gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

gl = axs[0,1].gridlines(draw_labels=False, crs=ccrs.PlateCarree(), linestyle='--', color='gray', alpha=0.7)
gl.top_labels = False  # Disable top labels
gl.right_labels = False  # Disable right labels
gl.bottom_labels=False
gl.xlabel_style = {'size': 10, 'color': 'black'}  # Customize longitude label style
gl.ylabel_style = {'size': 10, 'color': 'black'}  # Customize latitude label style

axs[0,0].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none'
)

axs[0,1].add_feature(
    cfeature.OCEAN.with_scale('50m'),
    facecolor='lightgray',
    edgecolor='none'
)

c1 = axs[0,0].pcolormesh(
    lon_sub, lat_sub, np.max(MLD_SO_3_sub, axis=0) - np.max(MLD_SO_1_sub, axis=0),
    transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-300, vmax=300)

cbar1 = fig.colorbar(c1, ax=axs[0, 0], orientation='vertical', extend='both', shrink=1)
cbar1.set_label('MLD maximum difference [m]', fontsize=12)
axs[0,0].set_title('a) MLD maximum difference', fontsize=14)

c2 = axs[0,1].pcolormesh(
    lon_sub, lat_sub, np.var(MLD_SO_3_sub, axis=0) - np.var(MLD_SO_1_sub, axis=0),transform=ccrs.PlateCarree(), cmap='seismic', shading='auto', vmin=-5000, vmax=5000)

cbar2 = fig.colorbar(c2, ax=axs[0,1], orientation='vertical', extend = 'both', shrink=1)
cbar2.set_label('MLD variance difference [m$^2$]', fontsize=12)
axs[0,1].set_title('b) MLD variance difference', fontsize=14)

#Define region for som regions
lon_min, lon_max = -50, 0
lat_min, lat_max = -50, -35

lon_min_NZ, lon_max_NZ = 160, 190
lat_min_NZ, lat_max_NZ = -70, -61

lon_min_AU, lon_max_AU = 80, 150
lat_min_AU, lat_max_AU = -70, -50

lon_min_PA, lon_max_PA = -110, -60
lat_min_PA, lat_max_PA = -70, -50

lon_min_p, lon_max_p = 170, 250
lat_min_p, lat_max_p = -60, -45

#WGKP
lon_min_wgkp, lon_max_wgkp = -35, 80
lat_min_wgkp, lat_max_wgkp = -80, -50

x_1	= np.arange(lon_min, lon_max + 0.1, 0.1)
y_1	= np.zeros(len(x_1)) + lat_min
y_2	= np.zeros(len(x_1)) + lat_max

y_3 = np.arange(lat_min, lat_max + 0.1, 0.1)
x_2 = np.zeros(len(y_3)) + lon_min
x_3 = np.zeros(len(y_3)) + lon_max

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

axs[0,1].plot(x_1_NZ, y_1_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_NZ, y_2_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_NZ, y_3_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_NZ, y_3_NZ, '-', color='magenta', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[0,1].plot(x_1_AU, y_1_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_AU, y_2_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_AU, y_3_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_AU, y_3_AU, '-', color='darkred', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

axs[0,1].plot(x_1_PA, y_1_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_PA, y_2_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_PA, y_3_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_PA, y_3_PA, '-', color='blue', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_p	= np.arange(lon_min_p, lon_max_p + 0.1, 0.1)
y_1_p	= np.zeros(len(x_1_p)) + lat_min_p
y_2_p	= np.zeros(len(x_1_p)) + lat_max_p

y_3_p = np.arange(lat_min_p, lat_max_p + 0.1, 0.1)
x_2_p = np.zeros(len(y_3_p)) + lon_min_p
x_3_p = np.zeros(len(y_3_p)) + lon_max_p

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp

axs[0,1].plot(x_1_wgkp, y_1_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_1_wgkp, y_2_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_2_wgkp, y_3_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)
axs[0,1].plot(x_3_wgkp, y_3_wgkp, '-g', linewidth = 5.0, transform=ccrs.PlateCarree(), zorder=10)


axs[1,0].plot(np.mean(PD_NZ_1, axis=0), depth, color='magenta', linewidth = 2.0, label='NZ')
axs[1,0].plot(np.mean(PD_WGKP_1, axis=0), depth, color='green', linewidth = 2.0, label='WGKP')
axs[1,0].plot(np.mean(PD_AU_1, axis=0), depth, color='darkred', linewidth = 2.0, label='AU')
axs[1,0].plot(np.mean(PD_PA_1, axis=0), depth, color='blue', linewidth = 2.0, label='PA')
axs[1,0].legend()
axs[1,0].set_ylim(2000, 1)  
axs[1,0].set_xlim(1026.4, 1027.8)
axs[1,0].set_title('c) Area-averaged PD (model year 1-100)', fontsize=14)
axs[1,0].set_xlabel('Potential density [kg/m$^3$]', fontsize=12)
axs[1,0].set_ylabel('Depth [m]', fontsize = 12)
axs[1,0].grid()

axs[1,1].plot(np.mean(PD_NZ_2, axis=0), depth, color='magenta', linewidth = 2.0, label='NZ')
axs[1,1].plot(np.mean(PD_WGKP_2, axis=0), depth, color='green', linewidth = 2.0, label='WGKP')
axs[1,1].plot(np.mean(PD_AU_2, axis=0), depth, color='darkred', linewidth = 2.0, label='AU')
axs[1,1].plot(np.mean(PD_PA_2, axis=0), depth, color='blue', linewidth = 2.0, label='PA')
axs[1,1].legend()
axs[1,1].set_ylim(2000, 1)  
axs[1,1].set_xlim(1026.4, 1027.8)
axs[1,1].set_title('d) Area-averaged PD (model year 500-600)', fontsize=14)
axs[1,1].set_xlabel('Potential density [kg/m$^3$]', fontsize=12)
axs[1,1].grid()

plt.tight_layout()
plt.savefig(directory_figures +'Figure_7_SOM_OS.pdf')
plt.show()


# %%
