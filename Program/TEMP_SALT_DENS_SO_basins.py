#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 10:01:47 2025

@author: 6008399

Density structures in different basins of the SO

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

    WGKP = len(depth)
    n0 = np.zeros(WGKP)

    for i in range(WGKP):
        if i == 0:
            # surface: (rho0 - rho1)/(z1 - z0)
            n0[i] = (density[i] - density[i+1]) / (depth[i+1] - depth[i])

        elif i == WGKP - 1:
            # bottom: (rho_{n-2} - rho_{n-1})/(z_{n-1}-z_{n-2})
            n0[i] = (density[i-1] - density[i]) / (depth[i] - depth[i-1])

        else:
            # centered: (rho[i-1] - rho[i+1]) / (z[i+1] - z[i-1])
            n0[i] = (density[i-1] - density[i+1]) / (depth[i+1] - depth[i-1])

    # Your formula: N² = -(g / rho0) * n0
    N2 = - (g / rho0) * n0
    
    if np.any(N2 < 0):
        print("N2 is unstable!")

    return n0, N2

# def compute_N2_from_profile(density, depth, rho0=1027.0, g=9.81):
#     density = np.asarray(density)
#     depth = np.asarray(depth)

#     drho_dz = np.gradient(density, depth)
#     N2 = (g / rho0) * drho_dz   # because depth increases downward

#     if np.any(N2 < 0):
#         print("N2 is unstable!")

#     return drho_dz, N2

def Welch(data_1, data_2):
    """Conducts Welch t-test; returns highest significance level achieved (0..1)."""

    mean_1 = np.mean(data_1)
    mean_2 = np.mean(data_2)

    std_1 = np.sqrt(1.0 / (len(data_1) - 1) * np.sum((data_1 - mean_1) ** 2.0))
    std_2 = np.sqrt(1.0 / (len(data_2) - 1) * np.sum((data_2 - mean_2) ** 2.0))

    t_welch = (mean_1 - mean_2) / np.sqrt((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2)))

    dof = (
        ((std_1**2.0 / len(data_1)) + (std_2**2.0 / len(data_2))) ** 2.0
        / (
            (std_1**4.0 / (len(data_1) ** 2.0 * (len(data_1) - 1)))
            + (std_2**4.0 / (len(data_2) ** 2.0 * (len(data_2) - 1)))
        )
    )

    sig_levels = np.arange(50, 100, 0.5) / 100.0
    t_crit = stats.t.ppf((1.0 + sig_levels) / 2.0, dof)

    sig_index = np.where(np.fabs(t_welch) > t_crit)[0]
    significant = 0.0
    if len(sig_index) > 0:
        significant = sig_levels[sig_index[-1]]

    return significant

#%% Mixed layer depth WGKP

fh = netcdf.Dataset(directory_zenodo + '/TEMP_SALT_DENS_year_1-100_zonal_averaged_-60E-25E_Atlantic_SO.nc','r')

time                = fh.variables['time'][:]   #time
depth               = fh.variables['depth'][:]  #depth
lat                 = fh.variables['lat'][:]    #latitude
temp_1_AT_time      = fh.variables['TEMP_all'][:]   #temperature
salt_1_AT_time      = fh.variables['SALT_all'][:]   #salinity
dens_1_AT_time      = fh.variables['PD_all'][:]     #potential density
temp_1_AT           = fh.variables['TEMP'][:]   #temperature
salt_1_AT           = fh.variables['SALT'][:]   #salinity
dens_1_AT           = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_zenodo + '/TEMP_SALT_DENS_year_500-600_zonal_averaged_-60E-25E_Atlantic_SO.nc','r')

time                = fh.variables['time'][:]   #time
depth               = fh.variables['depth'][:]  #depth
lat                 = fh.variables['lat'][:]    #latitude
temp_2_AT_time      = fh.variables['TEMP_all'][:]   #temperature
salt_2_AT_time      = fh.variables['SALT_all'][:]   #salinity
dens_2_AT_time      = fh.variables['PD_all'][:]     #potential density
temp_2_AT           = fh.variables['TEMP'][:]   #temperature
salt_2_AT           = fh.variables['SALT'][:]   #salinity
dens_2_AT           = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_1-100_zonal_averaged_60W_25E_Atlantic_transect_SO.nc','r')

depth          = fh.variables['depth'][:]  #depth
lat_u            = fh.variables['lat'][:]    #latitude
UVEL_1_AT      = fh.variables['U_VEL'][:]   #temperature

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_500-600_zonal_averaged_60W_25E_Atlantic_transect_SO.nc','r')

depth          = fh.variables['depth'][:]  #depth
lat_u            = fh.variables['lat'][:]    #latitude
UVEL_2_AT      = fh.variables['U_VEL'][:]   #temperature

fh.close()

fh = netcdf.Dataset(directory_zenodo + '/TEMP_SALT_DENS_year_1-100_zonal_averaged_lon_150-295_transect_SO_Pacific.nc','r')

time                = fh.variables['time'][:]   #time
depth               = fh.variables['depth'][:]  #depth
lat                 = fh.variables['lat'][:]    #latitude
temp_1_PA_time      = fh.variables['TEMP_all'][:]   #temperature
salt_1_PA_time      = fh.variables['SALT_all'][:]   #salinity
dens_1_PA_time      = fh.variables['PD_all'][:]     #potential density
temp_1_PA           = fh.variables['TEMP'][:]   #temperature
salt_1_PA           = fh.variables['SALT'][:]   #salinity
dens_1_PA           = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_zenodo + '/TEMP_SALT_DENS_year_500-600_zonal_averaged_lon_150-295_transect_SO_Pacific.nc','r')

time                = fh.variables['time'][:]   #time
depth               = fh.variables['depth'][:]  #depth
lat                 = fh.variables['lat'][:]    #latitude
temp_2_PA_time      = fh.variables['TEMP_all'][:]   #temperature
salt_2_PA_time      = fh.variables['SALT_all'][:]   #salinity
dens_2_PA_time      = fh.variables['PD_all'][:]     #potential density
temp_2_PA           = fh.variables['TEMP'][:]   #temperature
salt_2_PA           = fh.variables['SALT'][:]   #salinity
dens_2_PA           = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_1-100_zonal_averaged_lon_150-295_transect_SO_Pacific.nc','r')

depth          = fh.variables['depth'][:]  #depth
lat            = fh.variables['lat'][:]    #latitude
UVEL_1_PA      = fh.variables['U_VEL'][:]   #temperature

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_500-600_zonal_averaged_lon_150-295_transect_SO_Pacific.nc','r')

depth          = fh.variables['depth'][:]  #depth
lat            = fh.variables['lat'][:]    #latitude
UVEL_2_PA      = fh.variables['U_VEL'][:]   #temperature

fh.close()

fh = netcdf.Dataset(directory_zenodo + '/TEMP_SALT_DENS_year_1-100_zonal_averaged_25E-150E_Indian_SO.nc','r')

time                = fh.variables['time'][:]   #time
depth               = fh.variables['depth'][:]  #depth
lat                 = fh.variables['lat'][:]    #latitude
temp_1_IN_time      = fh.variables['TEMP_all'][:]   #temperature
salt_1_IN_time      = fh.variables['SALT_all'][:]   #salinity
dens_1_IN_time      = fh.variables['PD_all'][:]     #potential density
temp_1_IN           = fh.variables['TEMP'][:]   #temperature
salt_1_IN           = fh.variables['SALT'][:]   #salinity
dens_1_IN           = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_zenodo + '/TEMP_SALT_DENS_year_500-600_zonal_averaged_25E-150E_Indian_SO.nc','r')

time                = fh.variables['time'][:]   #time
depth               = fh.variables['depth'][:]  #depth
lat                 = fh.variables['lat'][:]    #latitude
temp_2_IN_time      = fh.variables['TEMP_all'][:]   #temperature
salt_2_IN_time      = fh.variables['SALT_all'][:]   #salinity
dens_2_IN_time      = fh.variables['PD_all'][:]     #potential density
temp_2_IN           = fh.variables['TEMP'][:]   #temperature
salt_2_IN           = fh.variables['SALT'][:]   #salinity
dens_2_IN           = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_1-100_zonal_averaged_25E_150E_Indian_transect_SO.nc','r')

depth          = fh.variables['depth'][:]  #depth
lat            = fh.variables['lat'][:]    #latitude
UVEL_1_IN      = fh.variables['U_VEL'][:]   #temperature

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_500-600_zonal_averaged_25E_150E_Indian_transect_SO.nc','r')

depth          = fh.variables['depth'][:]  #depth
lat            = fh.variables['lat'][:]    #latitude
UVEL_2_IN      = fh.variables['U_VEL'][:]   #temperature

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Volume_90S_10N_25E_150E_Indian_basin.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
lon_IN         = fh.variables['lon'][:]    #latitude
volume_IN      = fh.variables['volume'][:]   #salinity

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Volume_90S_10N_60W_25E_Atlantic_basin.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
lon_AT         = fh.variables['lon'][:]    #latitude
volume_AT      = fh.variables['volume'][:]   #salinity

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Volume_90S_10N_150E_295E_Pacific_basin.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
lon_PA         = fh.variables['lon'][:]    #latitude
volume_PA      = fh.variables['volume'][:]   #salinity

fh.close()

#%% Meridional density gradient

bin_size = 20
drho_dy1_AT = ma.masked_all((len(depth), len(lat)))
drho_dy2_AT = ma.masked_all((len(depth), len(lat)))
drho_dy1_PA = ma.masked_all((len(depth), len(lat)))
drho_dy2_PA = ma.masked_all((len(depth), len(lat)))
drho_dy1_IN = ma.masked_all((len(depth), len(lat)))
drho_dy2_IN = ma.masked_all((len(depth), len(lat)))

a1_AT = ma.masked_all((len(depth), len(lat)))

deg2m = 111e3   # degrees latitude → meters

for depth_i in range(len(depth)):
    for lat_i in range(bin_size, len(lat) - bin_size):

        lat_1 = lat_i - bin_size #southern
        lat_2 = lat_i + bin_size #northern
        
        print(lat[lat_1])
        print(lat[lat_2])
        
        #sys.exit()

        # actual latitude difference (in meters)
        del_lat_m = np.abs(lat[lat_2] - lat[lat_1])  * deg2m

        # density difference
        del_dens1_AT = dens_1_AT[depth_i, lat_2] - dens_1_AT[depth_i, lat_1]
        del_dens2_AT = dens_2_AT[depth_i, lat_2] - dens_2_AT[depth_i, lat_1]
        
        del_dens1_PA = dens_1_PA[depth_i, lat_2] - dens_1_PA[depth_i, lat_1]
        del_dens2_PA = dens_2_PA[depth_i, lat_2] - dens_2_PA[depth_i, lat_1]
        
        del_dens1_IN = dens_1_IN[depth_i, lat_2] - dens_1_IN[depth_i, lat_1]
        del_dens2_IN = dens_2_IN[depth_i, lat_2] - dens_2_IN[depth_i, lat_1]

        drho_dy1_AT[depth_i, lat_i] = del_dens1_AT / del_lat_m
        drho_dy2_AT[depth_i, lat_i] = del_dens2_AT / del_lat_m
        
        drho_dy1_PA[depth_i, lat_i] = del_dens1_PA / del_lat_m
        drho_dy2_PA[depth_i, lat_i] = del_dens2_PA / del_lat_m
        
        drho_dy1_IN[depth_i, lat_i] = del_dens1_IN / del_lat_m
        drho_dy2_IN[depth_i, lat_i] = del_dens2_IN / del_lat_m
        
        #Latitude range, southern boundary is set at y = 0 m
        del_lat_m = (lat[lat_1:lat_2+1] - lat[lat_1]) * deg2m
 
        #Density slope
        a, b  = np.polyfit(del_lat_m, dens_1_AT[depth_i, lat_1:lat_2+1], 1)
        a1_AT[depth_i, lat_i] = a
    
#%% Buoyancy frequency

n0_lat_1_AT = ma.masked_all((len(depth), len(lat)))
n0_lat_2_AT = ma.masked_all((len(depth), len(lat)))
N2_lat_1_AT = ma.masked_all((len(depth), len(lat)))
N2_lat_2_AT = ma.masked_all((len(depth), len(lat)))

n0_lat_1_PA = ma.masked_all((len(depth), len(lat)))
n0_lat_2_PA = ma.masked_all((len(depth), len(lat)))
N2_lat_1_PA = ma.masked_all((len(depth), len(lat)))
N2_lat_2_PA = ma.masked_all((len(depth), len(lat)))

n0_lat_1_IN = ma.masked_all((len(depth), len(lat)))
n0_lat_2_IN = ma.masked_all((len(depth), len(lat)))
N2_lat_1_IN = ma.masked_all((len(depth), len(lat)))
N2_lat_2_IN = ma.masked_all((len(depth), len(lat)))

for lat_i in range(len(lat)):
    n0_lat_1_AT[:, lat_i], N2_lat_1_AT[:,lat_i] = compute_N2_from_profile(dens_1_AT[:,lat_i], depth)
    n0_lat_2_AT[:, lat_i], N2_lat_2_AT[:,lat_i] = compute_N2_from_profile(dens_2_AT[:,lat_i], depth)
    
    n0_lat_1_PA[:, lat_i], N2_lat_1_PA[:,lat_i] = compute_N2_from_profile(dens_1_PA[:,lat_i], depth)
    n0_lat_2_PA[:, lat_i], N2_lat_2_PA[:,lat_i] = compute_N2_from_profile(dens_2_PA[:,lat_i], depth)
    
    n0_lat_1_IN[:, lat_i], N2_lat_1_IN[:,lat_i] = compute_N2_from_profile(dens_1_IN[:,lat_i], depth)
    n0_lat_2_IN[:, lat_i], N2_lat_2_IN[:,lat_i] = compute_N2_from_profile(dens_2_IN[:,lat_i], depth)
    

#%% Only until 30N

volume_lat_IN = np.sum(volume_IN[:,0:860,:], axis=2)
volume_lat_AT = np.sum(volume_AT[:,0:860,:], axis=2)
volume_lat_PA = np.sum(volume_PA[:,0:860,:], axis=2)

PD_weighted_IN_1 = np.sum(dens_1_IN[:,0:860] * volume_lat_IN, axis=1) / np.sum(volume_lat_IN, axis=1)
PD_weighted_IN_2 = np.sum(dens_2_IN[:,0:860] * volume_lat_IN, axis=1) / np.sum(volume_lat_IN, axis=1)

PD_weighted_AT_1 = np.sum(dens_1_AT[:,0:860] * volume_lat_AT, axis=1) / np.sum(volume_lat_AT, axis=1)
PD_weighted_AT_2 = np.sum(dens_2_AT[:,0:860] * volume_lat_AT, axis=1) / np.sum(volume_lat_AT, axis=1)

PD_weighted_PA_1 = np.sum(dens_1_PA[:,0:860] * volume_lat_PA, axis=1) / np.sum(volume_lat_PA, axis=1)
PD_weighted_PA_2 = np.sum(dens_2_PA[:,0:860] * volume_lat_PA, axis=1) / np.sum(volume_lat_PA, axis=1)

plt.figure()
plt.plot(PD_weighted_IN_1, depth, label='before')
plt.plot(PD_weighted_IN_2, depth, label='after')
plt.legend()
plt.ylim(2000, 0)

#%%

n0_som_1_IN, N2_som_1_IN = compute_N2_from_profile(PD_weighted_IN_1, depth)
n0_som_2_IN, N2_som_2_IN = compute_N2_from_profile(PD_weighted_IN_2, depth)

n0_som_1_AT, N2_som_1_AT = compute_N2_from_profile(PD_weighted_AT_1, depth)
n0_som_2_AT, N2_som_2_AT = compute_N2_from_profile(PD_weighted_AT_2, depth)

n0_som_1_PA, N2_som_1_PA = compute_N2_from_profile(PD_weighted_PA_1, depth)
n0_som_2_PA, N2_som_2_PA = compute_N2_from_profile(PD_weighted_PA_2, depth)

plt.figure()
plt.plot(N2_som_1_AT, depth)
plt.plot(N2_som_2_AT, depth)
plt.ylim(2000, 0)

plt.figure()
plt.plot(n0_som_1_AT, depth)
plt.plot(n0_som_2_AT, depth)
plt.ylim(2000, 0)


plt.figure()
plt.plot(N2_som_1_PA, depth)
plt.plot(N2_som_2_PA, depth)
plt.ylim(2000, 0)

plt.figure()
plt.plot(N2_som_1_IN, depth)
plt.plot(N2_som_2_IN, depth)
plt.ylim(2000, 0)

#%%

import matplotlib.pyplot as plt

# Create a 2x3 subplot
fig, axs = plt.subplots(2, 3, figsize=(10, 8), sharey=True)

# Common y-axis limits
ylim = (3000, 1)

# Top row: Absolute N^2
axs[0, 0].plot((N2_som_1_AT)/1e-6, depth, color='blue', label='(1-100)')
axs[0, 0].plot((N2_som_2_AT)/1e-6, depth, color='green', label='(500-600)')
axs[0, 0].set_title('a) Atlantic sector', fontsize=14)
axs[0, 0].set_xlabel('N$^2$ [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
axs[0, 0].set_ylabel('Depth [m]', fontsize=12)
axs[0, 0].legend(loc=4, fontsize=11)
axs[0, 0].set_ylim(*ylim)
axs[0, 0].set_xlim(-7, 67)
axs[0, 0].grid()

axs[0, 1].plot((N2_som_1_IN)/1e-6, depth, color='blue', label='(1-100)')
axs[0, 1].plot((N2_som_2_IN)/1e-6, depth, color='green', label='(500-600)')
axs[0, 1].set_title('b) Indian sector', fontsize=14)
axs[0, 1].set_xlabel('N$^2$ [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
axs[0, 1].set_ylim(*ylim)
axs[0, 1].set_xlim(-7, 67)
axs[0, 1].legend(loc=4, fontsize=11)
axs[0, 1].grid()

axs[0, 2].plot((N2_som_1_PA)/1e-6, depth, color='blue', label='(1-100)')
axs[0, 2].plot((N2_som_2_PA)/1e-6, depth, color='green', label='(500-600)')
axs[0, 2].set_title('c) Pacific sector', fontsize=14)
axs[0, 2].set_xlabel('N$^2$ [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
axs[0, 2].set_ylim(*ylim)
axs[0, 2].set_xlim(-7, 67)
axs[0, 2].legend(loc=4, fontsize=11)
axs[0, 2].grid()

# Bottom row: N^2 differences
axs[1, 0].plot((N2_som_2_AT - N2_som_1_AT)/1e-6, depth, color='red', label='(500-600) vs (1-100)')
axs[1, 0].vlines(x=0, ymin=ylim[0], ymax=ylim[1], color='black', linestyle='--')
axs[1, 0].set_title('d) Atlantic sector', fontsize=14)
axs[1, 0].set_xlabel('N$^2$ difference [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
axs[1, 0].set_ylabel('Depth [m]', fontsize=12)
axs[1, 0].set_ylim(*ylim)
axs[1, 0].set_xlim(-7, 7)
axs[1, 0].grid()

axs[1, 1].plot((N2_som_2_IN - N2_som_1_IN)/1e-6, depth, color='red', label='(500-600) vs (1-100)')
axs[1, 1].vlines(x=0, ymin=ylim[0], ymax=ylim[1], color='black', linestyle='--')
axs[1, 1].set_title('e) Indian sector', fontsize=14)
axs[1, 1].set_xlabel('N$^2$ difference [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
axs[1, 1].set_ylim(*ylim)
axs[1, 1].set_xlim(-7, 7)
axs[1, 1].grid()

axs[1, 2].plot((N2_som_2_PA - N2_som_1_PA)/1e-6, depth, color='red', label='(500-600) vs (1-100)')
axs[1, 2].vlines(x=0, ymin=ylim[0], ymax=ylim[1], color='black', linestyle='--')
axs[1, 2].set_title('f) Pacific sector', fontsize=14)
axs[1, 2].set_xlabel('N$^2$ difference [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
axs[1, 2].set_ylim(*ylim)
axs[1, 2].set_xlim(-7, 7)
axs[1, 2].grid()

# Adjust layout
plt.tight_layout()

# Save the figure
fig.savefig(directory_figures + 'N2_AT_PA_IN.pdf')

plt.show()

#%%

fig, axs = plt.subplots(1, 3, figsize=(10, 5))

ax1, ax2, ax3 = axs.flatten()

#ax1.plot(N2_som_2_WGKP - N2_som_1_WGKP, depth, color='blue', label='SOM2 vs SOM1')
ax1.plot((N2_som_2_AT - N2_som_1_AT)/1e-6, depth, color='red', label='(500-600) vs (1-100)')
ax1.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax1.set_title('a) Atlantic Sector', fontsize=14)
ax1.set_xlabel('N$^2$ difference [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
ax1.set_ylabel('Depth [m]', fontsize=12)
#ax1.legend(loc=4, fontsize=11)
ax1.set_ylim(5000, 1)
ax1.set_ylim(3000, 1)
ax1.set_xlim(-7, 7)
ax1.grid()

#ax2.plot(N2_som_2_WGKP - N2_som_1_WGKP, depth, color='blue', label='SOM2 vs SOM1')
ax2.plot((N2_som_2_IN - N2_som_1_IN)/1e-6, depth, color='red', label='(500-600) vs (1-100)')
ax2.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax2.set_title('b) Indian Sector', fontsize=14)
ax2.set_xlabel('N$^2$ difference [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
#ax2.set_ylabel('Depth [m]')
#ax2.legend(loc=4, fontsize=11)
ax2.set_ylim(5000, 1)
ax2.set_ylim(3000, 1)
ax2.set_xlim(-7., 7)
ax2.grid()

#ax2.plot(N2_som_2_WGKP - N2_som_1_WGKP, depth, color='blue', label='SOM2 vs SOM1')
ax3.plot((N2_som_2_PA - N2_som_1_PA)/1e-6, depth, color='red', label='(500-600) vs (1-100)')
ax3.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax3.set_title('c) Pacific Sector', fontsize=14)
ax3.set_xlabel('N$^2$ difference [x 1e$^{-6}$ s$^{-1}$]', fontsize=12)
#ax3.set_ylabel('Depth [m]')
#ax3.legend(loc=4, fontsize=11)
ax3.set_ylim(3000, 1)
#ax3.set_ylim(2000, 1)
ax3.set_xlim(-7., 7)
ax3.grid()

#plt.suptitle('Area-averaged N$^2$ difference SO30 (500-600) vs (1-100)', fontsize=16)
plt.tight_layout()

#fig.savefig(directory_figures +'N2_AT_PA_IN.pdf')

#%% Isopycnal slope

s1_AT =      - drho_dy1_AT/n0_lat_1_AT
s1_AT = np.ma.masked_where(N2_lat_1_AT < 1e-6, s1_AT)
s2_AT =     - drho_dy2_AT/n0_lat_2_AT
s2_AT = np.ma.masked_where(N2_lat_2_AT < 1e-6, s2_AT)

s1_PA =     - drho_dy1_PA/n0_lat_1_PA
s1_PA = np.ma.masked_where(N2_lat_1_PA < 1e-6, s1_PA)
s2_PA =    - drho_dy2_PA/n0_lat_2_PA
s2_PA = np.ma.masked_where(N2_lat_2_PA < 1e-6, s2_PA)

s1_IN =    - drho_dy1_IN/n0_lat_1_IN
s1_IN = np.ma.masked_where(N2_lat_1_IN < 1e-6, s1_IN)
s2_IN =     - drho_dy2_IN/n0_lat_2_IN
s2_IN = np.ma.masked_where(N2_lat_2_IN < 1e-6, s2_IN)
        
#%%

lat_idx = -1 #latitude of 0.05 degN (equator)
depth_idx0 = 18 #depth of 465m
depth_idx1 = 21 #depth of 918m
depth_idx2 = 25 #depth of 1380m
depth_idx3 = 30 #depth of 2125m

plt.figure(figsize=(12,6))
plt.title('Isopycnal slope')
plt.contourf(lat, depth, s2_AT - s1_AT, levels = np.arange(-0.01, 0.01, 0.001), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [1026.2], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [1026.2], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [1027.1], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [1027.1], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [1027.7], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [1027.7], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [1028.9], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [1028.9], linestyles = '--', colors='black')
plt.ylim(depth[-1], 0)
#plt.xlim(-80, -50)

#%%

plt.figure(figsize=(12,6))
plt.title('Density')
plt.contourf(lat, depth, dens_2_AT - dens_1_AT, levels = np.arange(-0.01, 0.01, 0.001), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1], 0)
#plt.xlim(-80, -50)

#%%
plt.figure(figsize=(12,6))
plt.title('a2 - a1')
plt.contourf(lat, depth, drho_dy2_AT - drho_dy1_AT, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1], 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('a1')
plt.contourf(lat, depth, drho_dy1_AT, levels = np.linspace(-0.000001, 0.000001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1], 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('a2')
plt.contourf(lat, depth, drho_dy2_AT, levels = np.linspace(-0.000001, 0.000001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1], 0)
#plt.xlim(-80, -50)

#%%

plt.figure(figsize=(12,6))
plt.title('s1')
plt.contourf(lat, depth, s1_AT, levels = np.linspace(-0.001, 0.001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('s2')
plt.contourf(lat, depth, s2_AT, levels = np.linspace(-0.001, 0.001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('s2 - s1')
plt.contourf(lat, depth, s2_AT - s1_AT, levels = np.linspace(-0.0005, 0.0005, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
#plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
#plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
#plt.contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
#plt.contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(2500, 0)
#plt.xlim(-50, -30)

#%%

# Create a 1x3 figure
fig, axs = plt.subplots(1, 3, figsize=(14, 4))  # 1 row, 3 columns

# First subplot: s1_AT
axs[0].set_title('a) Atlantic basin')
c1 = axs[0].contourf(lat, depth, s2_AT - s1_AT, levels=np.linspace(-0.001, 0.001, 21), extend='both', cmap='RdBu_r')
fig.colorbar(c1, ax=axs[0])
axs[0].contour(lat, depth, dens_1_AT, levels=[dens_1_AT[depth_idx0, lat_idx]], colors='black')
axs[0].contour(lat, depth, dens_2_AT, levels=[dens_2_AT[depth_idx0, lat_idx]], linestyles='--', colors='black')
axs[0].contour(lat, depth, dens_1_AT, levels=[dens_1_AT[depth_idx1, lat_idx]], colors='black')
axs[0].contour(lat, depth, dens_2_AT, levels=[dens_2_AT[depth_idx1, lat_idx]], linestyles='--', colors='black')
axs[0].set_ylim(3000, 0)
axs[0].set_ylabel('Depth [m]', fontsize=12)
axs[0].set_xlabel('Latitude', fontsize=12)

# Second subplot: s2_AT
axs[1].set_title('b) Indian basin')
c2 = axs[1].contourf(lat, depth, s2_IN - s1_IN, levels=np.linspace(-0.001, 0.001, 21), extend='both', cmap='RdBu_r')
fig.colorbar(c2, ax=axs[1])
axs[1].contour(lat, depth, dens_1_IN, levels=[dens_1_IN[depth_idx0, lat_idx]], colors='black')
axs[1].contour(lat, depth, dens_2_IN, levels=[dens_2_IN[depth_idx0, lat_idx]], linestyles='--', colors='black')
axs[1].contour(lat, depth, dens_1_IN, levels=[dens_1_IN[depth_idx1, lat_idx]], colors='black')
axs[1].contour(lat, depth, dens_2_IN, levels=[dens_2_IN[depth_idx1, lat_idx]], linestyles='--', colors='black')
axs[1].set_ylim(3000, 0)
axs[1].set_xlabel('Latitude', fontsize=12)

# Third subplot: s2_AT - s1_AT
axs[2].set_title('c) Pacific basin')
c3 = axs[2].contourf(lat, depth, s2_PA - s1_PA, levels=np.linspace(-0.0005, 0.0005, 21), extend='both', cmap='RdBu_r')
fig.colorbar(c3, ax=axs[2])
axs[2].contour(lat, depth, dens_1_PA, levels=[dens_1_PA[depth_idx0, lat_idx]], colors='black')
axs[2].contour(lat, depth, dens_2_PA, levels=[dens_2_PA[depth_idx0, lat_idx]], linestyles='--', colors='black')
axs[2].contour(lat, depth, dens_1_PA, levels=[dens_1_PA[depth_idx1, lat_idx]], colors='black')
axs[2].contour(lat, depth, dens_2_PA, levels=[dens_2_PA[depth_idx1, lat_idx]], linestyles='--', colors='black')
axs[2].set_ylim(3000, 0)
axs[2].set_xlabel('Latitude', fontsize=12)

# Adjust layout
plt.suptitle('Difference isopycnal slope (model year 500--600 minus 1--100)')
plt.tight_layout()
plt.savefig(directory_figures +'diff_isopycnal_slope_SO30_sectors.pdf')
plt.show()

#%%

plt.figure(figsize=(12,6))
plt.title('a')
plt.contourf(lat, depth, s1_PA, levels = np.linspace(-0.001, 0.001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('a')
plt.contourf(lat, depth, s2_PA, levels = np.linspace(-0.001, 0.001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('Isopycnal slope')
plt.contourf(lat, depth, s2_PA - s1_PA, levels = np.linspace(-0.0005, 0.0005, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

#%%

plt.figure(figsize=(12,6))
plt.title('a')
plt.contourf(lat, depth, s1_IN, levels = np.linspace(-0.001, 0.001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('a')
plt.contourf(lat, depth, s2_IN, levels = np.linspace(-0.001, 0.001, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

plt.figure(figsize=(12,6))
plt.title('Isopycnal slope')
plt.contourf(lat, depth, s2_IN - s1_IN, levels = np.linspace(-0.0005, 0.0005, 21), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(2500, 0)
#plt.xlim(-80, -50)

#%%

lat_idx = -1 #latitude of 0.05 degN (equator)
depth_idx0 = 10 #depth of 465m
depth_idx1 = 20 #depth of 918m
depth_idx2 = 25 #depth of 1380m
depth_idx3 = 30 #depth of 2125m

plt.figure()
plt.contourf(lat, depth, temp_1_AT)

plt.figure()
plt.contourf(lat, depth, salt_1_AT)
plt.colorbar()
plt.contour(lat, depth, salt_1_AT, levels = [salt_1_AT[depth_idx0, lat_idx]], colors = 'black')
#plt.contour(lat, depth, salt_2_AT, levels = [salt_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, salt_1_AT, levels = [salt_1_AT[depth_idx1, lat_idx]], colors = 'black')
#axs[0,1].contour(lat, depth, salt_2_AT, levels = [salt_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, salt_1_AT, levels = [salt_1_AT[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, salt_1_AT, levels = [salt_1_AT[depth_idx3, lat_idx]], colors = 'black')
#axs[0,1].contour(lat, depth, salt_2_AT, levels = [salt_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1],0)

#%%
plt.figure()
plt.contourf(lat, depth, UVEL_1_PA, levels = np.linspace(-0.4, 0.4, 21), cmap='seismic')
plt.colorbar()
#plt.contour(lat_u, depth, UVEL_1_AT, levels= [0.05])#np.linspace(0.03, 0.05,2))
#plt.colorbar()
plt.ylim(depth[-1],0)

#%% Atlantic sector


divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=6)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2_AT - temp_1_AT, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')

#axs[0,0].contour(lat, depth, temp_1_AT, levels = [temp_1_AT[depth_idx0, lat_idx]], colors = 'black')
#axs[0,0].contour(lat, depth, temp_2_AT, levels = [temp_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
#axs[0,0].contour(lat, depth, temp_1_AT, levels = [temp_1_AT[depth_idx1, lat_idx]], colors = 'black')
#axs[0,0].contour(lat, depth, temp_2_AT, levels = [temp_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
#axs[0,0].contour(lat, depth, temp_1_AT, levels = [temp_1_AT[depth_idx3, lat_idx]], colors = 'black')
#axs[0,0].contour(lat, depth, temp_2_AT, levels = [temp_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')

axs[0,0].set_xlim(-75,8)
axs[0,0].set_ylim(2500, 0)
#axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,0].tick_params(labelsize=11)
axs[0,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,0].set_title('a) Temperature', fontsize=15)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(temp_1_AT_time[:, depth_i, lat_i], temp_2_AT_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,0].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=6, facecolors="none")

CS2 = axs[0,1].contourf(lat, depth, salt_2_AT - salt_1_AT, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,8)
axs[0,1].set_ylim(2500, 0)
#axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,1].tick_params(labelsize=11)
axs[0,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,1].set_title('b) Salinity', fontsize=15)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(salt_1_AT_time[:, depth_i, lat_i], salt_2_AT_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,1].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=6, facecolors="none")


CS3 = axs[1,0].contourf(lat, depth,  (drho_dy2_AT - drho_dy1_AT)/1e-7, levels = np.linspace(-0.0000005/1e-7, 0.0000005/1e-7, 21), extend = 'both', cmap = 'PuOr_r')
# contour1 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[dens_1_AT[depth_idx0, lat_idx]], colors='black')
# contour2 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[dens_2_AT[depth_idx0, lat_idx]], linestyles='--', colors='black')
# contour3 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[dens_1_AT[depth_idx1, lat_idx]], colors='black')
# contour4 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[dens_2_AT[depth_idx1, lat_idx]], linestyles='--', colors='black')
# contour5 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[dens_1_AT[depth_idx2, lat_idx]], colors='black')
# contour6 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[dens_2_AT[depth_idx2, lat_idx]], linestyles='--', colors='black')

contour1 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[1027.0], colors='black')
contour2 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[1027.5], colors='black')
contour4 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  0].contour(lat, depth, dens_1_AT, levels=[1027.7], colors='black')
contour6 = axs[1,  0].contour(lat, depth, dens_2_AT, levels=[1027.7], linestyles='--', colors='black')

axs[1, 0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005/1e-7, 0.00000051/1e-7, 0.0000005/1e-7))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta \rho / \Delta y$ difference [x 10$^-7$ kg/m$^4$]', fontsize = 12)

axs[1,0].set_title('c) Meridional density gradient', fontsize=15)
axs[1,0].set_xlim(-75,8)
axs[1,0].set_ylim(2500, 0)
#axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,0].tick_params(labelsize=11)
axs[1,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

CS2 = axs[1,1].contourf(lat_u, depth, UVEL_2_AT - UVEL_1_AT, levels= np.linspace(-0.05, 0.05, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.05, 0.051, 0.05))
cbar.set_label(r'UVEL difference [m/s]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[1,1].set_title('d) Zonal velocity', fontsize=15)
axs[1,1].set_xlim(-75, 8)
axs[1,1].set_ylim(2500, 0)
#axs[1,1].contour(lat_u, depth, UVEL_1_AT, levels= [0.05], colors = 'black')
axs[1, 1].axvline(x=-55, linestyle='--', color='black')
axs[1, 1].axvline(x=-35, linestyle='--', color='black')
#axs[1,1].contour(lat_u, depth, UVEL_2_AT, levels= [-0.03, 0.03], colors = 'grey')
axs[1,1].tick_params(labelsize=11)
axs[1,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
#axs[1,1].set_ylabel('Depth [m]', fontsize=13)
#axs[1,1].tick_params(axis='x', labelsize=11)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

#plt.suptitle('Atlantic sector', fontsize=19)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_3_SOM_OS.pdf')
plt.show()

#%% See if density anomalies are more related to temperature or salinity contributions. Linearly decompose density into its salinity and temperature contributions

def RHO_0_dT(T, S):
    """Temperature driven density changes. EOS derived from (Millero and Poisson, 1981)"""
    #Reference density which is not pressure dependent
    a_0 = 6.793952 * 10**(-2.0) 
    a_1 = -9.095290 * 10**(-3.0) * 2
    a_2 = 1.001685 * 10**(-4.0) * 3 
    a_3 = - 1.120083 * 10**(-6.0) * 4
    a_4 = 6.536332 * 10**(-9.0) * 5
 
    b_0 = - 4.4490 * 10**(-3.0)
    b_1 = 1.0485 * 10**(-4.0) * 2
    b_2 =  - 1.2580 * 10**(-6.0) * 3
    b_3 =  3.315 * 10**(-9.0) * 4
 
    c_0 =  2.8441 * 10**(-4.0)
    c_1 = - 1.6871 * 10**(-5.0) * 2
    c_2 =  2.83258 * 10**(-7.0) * 3
 
    d_0 = - 1.97975 * 10**(-5.0)
    d_1 = 1.6641 * 10**(-6.0) * 2
    d_2 = - 3.1203 * 10**(-8.0) * 3
    rho_dT = (a_0+a_1*T+a_2*T**2+a_3*T**3+a_4*T**4)+(b_0+b_1*T+b_2*T**2+b_3*T**3)*S+(c_0+c_1*T+c_2*T**2)*S**(3/2)+(d_0+d_1*T+d_2*T**2)*S**2
    return rho_dT

def RHO_0_dS(T, S):
    """Salinity driven density changes. EOS derived from (Millero and Poisson, 1981)"""
    a_0 = 8.25917 * 10**(-1.0)
    a_1 = - 6.33761 * 10**(-3.0)
    a_2 = 5.4705 * 10**(-4.0)
    b_0 = - 4.4490 * 10**(-3.0)
    b_1 = 1.0485 * 10**(-4.0) 
    b_2 =  - 1.2580 * 10**(-6.0) 
    b_3 =  3.315 * 10**(-9.0)
 
    c_0 =  2.8441 * 10**(-4.0)
    c_1 = - 1.6871 * 10**(-5.0) 
    c_2 =  2.83258 * 10**(-7.0)
 
    d_0 = - 1.97975 * 10**(-5.0)
    d_1 = 1.6641 * 10**(-6.0) 
    d_2 = - 3.1203 * 10**(-8.0)
 
    rho_dS = (a_0+b_0*T+b_1*T**2+b_2*T**3+b_3*T**4)+(3/2)*(a_1+c_0*T+c_1*T**2+c_2*T**3)*S**(1/2)+2*(a_2+d_0*T+d_1*T**2+d_2*T**3)*S
    return rho_dS

#%%

#Differences
dT_AT = temp_2_AT - temp_1_AT
dS_AT = salt_2_AT - salt_1_AT
dR_AT = dens_2_AT - dens_1_AT

dT_PA = temp_2_PA - temp_1_PA
dS_PA = salt_2_PA - salt_1_PA
dR_PA = dens_2_PA - dens_1_PA

dT_IN = temp_2_IN - temp_1_IN
dS_IN = salt_2_IN - salt_1_IN
dR_IN = dens_2_IN - dens_1_IN

#Reference state for the derivatives (midpoint state)
T_ref_AT = 0.5 * (temp_1_AT + temp_2_AT)
S_ref_AT = 0.5 * (salt_1_AT + salt_2_AT)

T_ref_PA = 0.5 * (temp_1_PA + temp_2_PA)
S_ref_PA = 0.5 * (salt_1_PA + salt_2_PA)

T_ref_IN = 0.5 * (temp_1_IN + temp_2_IN)
S_ref_IN = 0.5 * (salt_1_IN + salt_2_IN)

drho_dT_AT = RHO_0_dT(T_ref_AT, S_ref_AT)   # shape: (depth, lat)
drho_dS_AT = RHO_0_dS(T_ref_AT, S_ref_AT)

drho_dT_PA = RHO_0_dT(T_ref_PA, S_ref_PA)   # shape: (depth, lat)
drho_dS_PA = RHO_0_dS(T_ref_PA, S_ref_PA)

drho_dT_IN = RHO_0_dT(T_ref_IN, S_ref_IN)   # shape: (depth, lat)
drho_dS_IN = RHO_0_dS(T_ref_IN, S_ref_IN)

#Decompose density change
dR_T_AT = drho_dT_AT * dT_AT    # temperature-driven density change
dR_S_AT = drho_dS_AT * dS_AT    # salinity-driven density change

dR_T_PA = drho_dT_PA * dT_PA    # temperature-driven density change
dR_S_PA = drho_dS_PA * dS_PA    # salinity-driven density change

dR_T_IN = drho_dT_IN * dT_IN    # temperature-driven density change
dR_S_IN = drho_dS_IN * dS_IN    # salinity-driven density change

#linear reconstructed total density change
dR_lin_AT = dR_T_AT + dR_S_AT
dR_lin_PA = dR_T_PA + dR_S_PA
dR_lin_IN = dR_T_IN + dR_S_IN

#Reconstruction error
recon_error_AT = dR_AT - dR_lin_AT
recon_error_PA = dR_PA - dR_lin_PA
recon_error_IN = dR_IN - dR_lin_IN

# %%

g = 9.81
rho0 = 1027.0

def d_dz(var, depth):
    return np.gradient(var, depth, axis=0)

# Actual N² in each state for depth increasing downward!! (so positive N² means stable stratification)
N2_1_AT = (g/rho0) * d_dz(dens_1_AT, depth)
N2_2_AT = (g/rho0) * d_dz(dens_2_AT, depth)
dN2_AT = N2_lat_2_AT - N2_lat_1_AT

N2_1_IN = (g/rho0) * d_dz(dens_1_IN, depth)
N2_2_IN = (g/rho0) * d_dz(dens_2_IN, depth)
dN2_IN = N2_lat_2_IN - N2_lat_1_IN

N2_1_PA = (g/rho0) * d_dz(dens_1_PA, depth)
N2_2_PA = (g/rho0) * d_dz(dens_2_PA, depth)
dN2_PA = N2_lat_2_PA - N2_lat_1_PA

# T/S contributions to ΔN²
dN2_T_AT = (g/rho0) * d_dz(dR_T_AT, depth)
dN2_S_AT = (g/rho0) * d_dz(dR_S_AT, depth)
dN2_lin_AT = dN2_T_AT + dN2_S_AT

dN2_T_IN = (g/rho0) * d_dz(dR_T_IN, depth)
dN2_S_IN = (g/rho0) * d_dz(dR_S_IN, depth)
dN2_lin_IN = dN2_T_IN + dN2_S_IN

dN2_T_PA = (g/rho0) * d_dz(dR_T_PA, depth)
dN2_S_PA = (g/rho0) * d_dz(dR_S_PA, depth)
dN2_lin_PA = dN2_T_PA + dN2_S_PA    

# Reconstruction error
recon_error_N2_AT = dN2_AT - dN2_lin_AT
recon_error_N2_IN = dN2_IN - dN2_lin_IN
recon_error_N2_PA = dN2_PA - dN2_lin_PA

#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=6)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(3, 2, figsize=(14, 10))

#Temperature
CS = axs[0,0].contourf(lat, depth, temp_2_AT - temp_1_AT, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')

axs[0,0].set_xlim(-75,8)
axs[0,0].set_ylim(2500, 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,0].tick_params(labelsize=12)
#axs[0,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[0,0].set_ylabel('Depth [m]', fontsize=15)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,0].set_title('a) Temperature', fontsize=17)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(temp_1_AT_time[:, depth_i, lat_i], temp_2_AT_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,0].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=15, facecolors="none")

#Salinity
CS2 = axs[0,1].contourf(lat, depth, salt_2_AT - salt_1_AT, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,8)
axs[0,1].set_ylim(2500, 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,1].tick_params(labelsize=12)
#axs[0,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,1].set_title('b) Salinity', fontsize=17)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(salt_1_AT_time[:, depth_i, lat_i], salt_2_AT_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,1].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=15, facecolors="none")


#Density
CS3 = axs[1,0].contourf(lat, depth,  dens_2_AT - dens_1_AT, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

contour1 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[1027.0], colors='black')
contour2 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1, 0].contour(lat, depth, dens_1_AT, levels=[1027.5], colors='black')
contour4 = axs[1, 0].contour(lat, depth, dens_2_AT, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  0].contour(lat, depth, dens_1_AT, levels=[1027.7], colors='black')
contour6 = axs[1,  0].contour(lat, depth, dens_2_AT, levels=[1027.7], linestyles='--', colors='black')

axs[1, 0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.5, 0.51, 0.2))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'Density difference [kg/m$^3$]', fontsize = 12)

axs[1,0].set_title('c) Potential density', fontsize=17)
axs[1,0].set_xlim(-75,8)
axs[1,0].set_ylim(2500, 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,0].tick_params(labelsize=12)
#axs[1,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,0].set_ylabel('Depth [m]', fontsize=15)
#axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

#Meridional density gradient
CS3 = axs[1,1].contourf(lat, depth,  (drho_dy2_AT - drho_dy1_AT)/1e-7, levels = np.linspace(-0.0000005/1e-7, 0.0000005/1e-7, 21), extend = 'both', cmap = 'PuOr_r')

contour1 = axs[1, 1].contour(lat, depth, dens_1_AT, levels=[1027.0], colors='black')
contour2 = axs[1, 1].contour(lat, depth, dens_2_AT, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1, 1].contour(lat, depth, dens_1_AT, levels=[1027.5], colors='black')
contour4 = axs[1, 1].contour(lat, depth, dens_2_AT, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  1].contour(lat, depth, dens_1_AT, levels=[1027.7], colors='black')
contour6 = axs[1,  1].contour(lat, depth, dens_2_AT, levels=[1027.7], linestyles='--', colors='black')

axs[1, 1].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005/1e-7, 0.00000051/1e-7, 0.0000005/1e-7))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta \rho / \Delta y$ difference [x 10$^-7$ kg/m$^4$]', fontsize = 12)

axs[1,1].set_title('d) Meridional density gradient', fontsize=17)
axs[1,1].set_xlim(-75,8)
axs[1,1].set_ylim(2500, 0)
axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,1].tick_params(labelsize=12)
#axs[1,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,1].set_ylabel('Depth [m]', fontsize=15)
#axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=14)


#Stratification frequency
CS3 = axs[2,0].contourf(lat, depth,  N2_lat_2_AT - N2_lat_1_AT, levels= np.linspace(-5e-6, 5e-6, 21), extend = 'both', cmap = 'RdBu_r')

contour1 = axs[2, 0].contour(lat, depth, dens_1_AT, levels=[1027.0], colors='black')
contour2 = axs[2, 0].contour(lat, depth, dens_2_AT, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[2, 0].contour(lat, depth, dens_1_AT, levels=[1027.5], colors='black')
contour4 = axs[2, 0].contour(lat, depth, dens_2_AT, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[2,  0].contour(lat, depth, dens_1_AT, levels=[1027.7], colors='black')
contour6 = axs[2,  0].contour(lat, depth, dens_2_AT, levels=[1027.7], linestyles='--', colors='black')

axs[2, 0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-5e-6, 5e-6 + 1e-7, 1e-6))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 12)

axs[2,0].set_title('e) N$^2$', fontsize=17)
axs[2,0].set_xlim(-75,8)
axs[2,0].set_ylim(2500, 0)
axs[2,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[2,0].tick_params(labelsize=12)
#axs[2,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[2,0].set_ylabel('Depth [m]', fontsize=15)
axs[2,0].set_xlabel('Latitude [$^\circ$S]', fontsize=15)

#Stratification frequency
CS3 = axs[2,1].contourf(lat, depth,  dN2_S_AT, levels= np.linspace(-5e-6, 5e-6, 21), extend = 'both', cmap = 'RdBu_r')

contour1 = axs[2, 1].contour(lat, depth, dens_1_AT, levels=[1027.0], colors='black')
contour2 = axs[2, 1].contour(lat, depth, dens_2_AT, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[2, 1].contour(lat, depth, dens_1_AT, levels=[1027.5], colors='black')
contour4 = axs[2, 1].contour(lat, depth, dens_2_AT, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[2,  1].contour(lat, depth, dens_1_AT, levels=[1027.7], colors='black')
contour6 = axs[2,  1].contour(lat, depth, dens_2_AT, levels=[1027.7], linestyles='--', colors='black')

axs[2, 1].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-5e-6, 5e-6 + 1e-7, 1e-6))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 12)

axs[2,1].set_title('f) Salinity-driven N$^2$', fontsize=17)
axs[2,1].set_xlim(-75,8)
axs[2,1].set_ylim(2500, 0)
axs[2,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[2,1].tick_params(labelsize=12)
#axs[2,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[2,1].set_ylabel('Depth [m]', fontsize=15)
axs[2,1].set_xlabel('Latitude [$^\circ$S]', fontsize=15)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_3_SOM_OS_revised.pdf')
plt.show()

#%% Indian sector

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2_IN - temp_1_IN, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')

axs[0,0].set_xlim(-70,8)
axs[0,0].set_ylim(2500, 0)
axs[0,0].set_xticklabels(['70', '60', '50', '40',  '30', '20', '10', '0', '-10'])
axs[0,0].tick_params(labelsize=11)
#axs[0,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,0].set_title('a) Temperature', fontsize=15)

CS2 = axs[0,1].contourf(lat, depth, salt_2_IN - salt_1_IN, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-70,8)
axs[0,1].set_ylim(2500, 0)
axs[0,1].set_xticklabels(['70', '60', '50', '40',  '30', '20', '10', '0', '-10'])
axs[0,1].tick_params(labelsize=11)
#axs[0,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,1].set_title('b) Salinity', fontsize=15)


CS3 = axs[1,0].contourf(lat, depth,  (drho_dy2_IN - drho_dy1_IN)/1e-7, levels = np.linspace(-0.0000005/1e-7, 0.0000005/1e-7, 21), extend = 'both', cmap = 'PuOr_r')
contour1 = axs[1,0].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx0, lat_idx]], colors = 'black')
contour2 = axs[1,0].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx0, lat_idx]], linestyles = '--', colors='black')
contour3 = axs[1,0].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx1, lat_idx]], colors = 'black')
contour4 = axs[1,0].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx1, lat_idx]], linestyles = '--', colors='black')
contour5 = axs[1,0].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx2, lat_idx]], colors = 'black')
contour6 = axs[1,0].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx2, lat_idx]], linestyles = '--', colors='black')
#plt.legend()

cbar	= colorbar(CS3, ticks = np.arange(-0.0000005/1e-7, 0.00000051/1e-7, 0.0000005/1e-7))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta \rho / \Delta y$ difference [x 10$^-7$ kg/m$^4$]', fontsize = 12)

axs[1,0].set_title('c) Meridional density gradient', fontsize=15)
axs[1,0].set_xlim(-70,8)
axs[1,0].set_ylim(2500, 0)
axs[1,0].set_xticklabels(['70', '60', '50', '40',  '30', '20', '10', '0', '-10'])
axs[1,0].tick_params(labelsize=11)
#axs[1,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

CS2 = axs[1,1].contourf(lat, depth, UVEL_2_IN - UVEL_1_IN, levels= np.linspace(-0.05, 0.05, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.05, 0.051, 0.05))
cbar.set_label(r'UVEL difference [m/s]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[1,1].set_title('d) Zonal velocity', fontsize=15)
axs[1,1].set_xlim(-70, 8)
axs[1,1].set_ylim(2500, 0)
#axs[1,1].contour(lat, depth, UVEL_1_IN, levels= [0.05], colors = 'black')
axs[1, 1].axvline(x=-55, linestyle='--', color='black')
axs[1, 1].axvline(x=-35, linestyle='--', color='black')
axs[1,1].tick_params(labelsize=11)
axs[1,1].set_xticklabels(['70', '60', '50', '40',  '30', '20', '10', '0', '-10'])
#axs[1,1].set_ylabel('Depth [m]', fontsize=13)
#axs[1,1].tick_params(axis='x', labelsize=11)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

plt.suptitle('Indian sector', fontsize=19)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_S2_SOM_OS.pdf')
plt.show()

#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=6)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(3, 2, figsize=(14, 10))

#Temperature
CS = axs[0,0].contourf(lat, depth, temp_2_IN - temp_1_IN, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')

axs[0,0].set_xlim(-70,8)
axs[0,0].set_ylim(2500, 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,0].tick_params(labelsize=12)
#axs[0,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[0,0].set_ylabel('Depth [m]', fontsize=15)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,0].set_title('a) Temperature', fontsize=17)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(temp_1_IN_time[:, depth_i, lat_i], temp_2_IN_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,0].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=15, facecolors="none")

#Salinity
CS2 = axs[0,1].contourf(lat, depth, salt_2_IN - salt_1_IN, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-70,8)
axs[0,1].set_ylim(2500, 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,1].tick_params(labelsize=12)
#axs[0,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,1].set_title('b) Salinity', fontsize=17)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(salt_1_IN_time[:, depth_i, lat_i], salt_2_IN_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,1].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=15, facecolors="none")


#Density
CS3 = axs[1,0].contourf(lat, depth,  dens_2_IN - dens_1_IN, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

contour1 = axs[1, 0].contour(lat, depth, dens_1_IN, levels=[1027.0], colors='black')
contour2 = axs[1, 0].contour(lat, depth, dens_2_IN, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1,  0].contour(lat, depth, dens_1_IN, levels=[1027.5], colors='black')
contour4 = axs[1,  0].contour(lat, depth, dens_2_IN, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  0].contour(lat, depth, dens_1_IN, levels=[1027.7], colors='black')
contour6 = axs[1,  0].contour(lat, depth, dens_2_IN, levels=[1027.7], linestyles='--', colors='black')

axs[1,  0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.5, 0.51, 0.2))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'Density difference [kg/m$^3$]', fontsize = 12)

axs[1,0].set_title('c) Potential density', fontsize=17)
axs[1,0].set_xlim(-70,8)
axs[1,0].set_ylim(2500, 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,0].tick_params(labelsize=12)
#axs[1,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,0].set_ylabel('Depth [m]', fontsize=15)
#axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

#Meridional density gradient
CS3 = axs[1,1].contourf(lat, depth,  (drho_dy2_IN - drho_dy1_IN)/1e-7, levels = np.linspace(-0.0000005/1e-7, 0.0000005/1e-7, 21), extend = 'both', cmap = 'PuOr_r')

contour1 = axs[1, 1].contour(lat, depth, dens_1_IN, levels=[1027.0], colors='black')
contour2 = axs[1, 1].contour(lat, depth, dens_2_IN, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1, 1].contour(lat, depth, dens_1_IN, levels=[1027.5], colors='black')
contour4 = axs[1,  1].contour(lat, depth, dens_2_IN, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  1].contour(lat, depth, dens_1_IN, levels=[1027.7], colors='black')
contour6 = axs[1,  1].contour(lat, depth, dens_2_IN, levels=[1027.7], linestyles='--', colors='black')

axs[1, 1].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  1].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  1].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  1].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  1].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1,  1].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005/1e-7, 0.00000051/1e-7, 0.0000005/1e-7))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta \rho / \Delta y$ difference [x 10$^-7$ kg/m$^4$]', fontsize = 12)

axs[1,1].set_title('d) Meridional density gradient', fontsize=17)
axs[1,1].set_xlim(-70,8)
axs[1,1].set_ylim(2500, 0)
axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,1].tick_params(labelsize=12)
#axs[1,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,1].set_ylabel('Depth [m]', fontsize=15)
#axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=14)


#Stratification frequency
CS3 = axs[2,0].contourf(lat, depth,  N2_lat_2_IN - N2_lat_1_IN, levels= np.linspace(-5e-6, 5e-6, 21), extend = 'both', cmap = 'RdBu_r')

contour1 = axs[2, 0].contour(lat, depth, dens_1_IN, levels=[1027.0], colors='black')
contour2 = axs[2, 0].contour(lat, depth, dens_2_IN, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[2,  0].contour(lat, depth, dens_1_IN, levels=[1027.5], colors='black')
contour4 = axs[2,  0].contour(lat, depth, dens_2_IN, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[2,  0].contour(lat, depth, dens_1_IN, levels=[1027.7], colors='black')
contour6 = axs[2,  0].contour(lat, depth, dens_2_IN, levels=[1027.7], linestyles='--', colors='black')

axs[2,  0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-5e-6, 5e-6 + 1e-7, 1e-6))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 12)

axs[2,0].set_title('e) N$^2$', fontsize=17)
axs[2,0].set_xlim(-70,8)
axs[2,0].set_ylim(2500, 0)
axs[2,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[2,0].tick_params(labelsize=12)
#axs[2,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[2,0].set_ylabel('Depth [m]', fontsize=15)
axs[2,0].set_xlabel('Latitude [$^\circ$S]', fontsize=15)

#Stratification frequency
CS3 = axs[2,1].contourf(lat, depth,  dN2_S_IN, levels= np.linspace(-5e-6, 5e-6, 21), extend = 'both', cmap = 'RdBu_r')

contour1 = axs[2, 1].contour(lat, depth, dens_1_IN, levels=[1027.0], colors='black')
contour2 = axs[2, 1].contour(lat, depth, dens_2_IN, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[2, 1].contour(lat, depth, dens_1_IN, levels=[1027.5], colors='black')
contour4 = axs[2,  1].contour(lat, depth, dens_2_IN, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[2,  1].contour(lat, depth, dens_1_IN, levels=[1027.7], colors='black')
contour6 = axs[2,  1].contour(lat, depth, dens_2_IN, levels=[1027.7], linestyles='--', colors='black')

axs[2,  1].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  1].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  1].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  1].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  1].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[2,  1].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-5e-6, 5e-6 + 1e-7, 1e-6))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 12)

axs[2,1].set_title('f) Salinity-driven N$^2$', fontsize=17)
axs[2,1].set_xlim(-70,8)
axs[2,1].set_ylim(2500, 0)
axs[2,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[2,1].tick_params(labelsize=12)
#axs[2,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[2,1].set_ylabel('Depth [m]', fontsize=15)
axs[2,1].set_xlabel('Latitude [$^\circ$S]', fontsize=15)

plt.suptitle('Indian sector', fontsize=19)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_S2_SOM_OS_revised.pdf')
plt.show()

#%% Pacific sector

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2_PA - temp_1_PA, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')

axs[0,0].set_xlim(-75,8)
axs[0,0].set_ylim(2500, 0)
#axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,0].tick_params(labelsize=11)
axs[0,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[0,0].set_ylabel('Depth [m]', fontsize=14)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,0].set_title('a) Temperature', fontsize=15)

CS2 = axs[0,1].contourf(lat, depth, salt_2_PA - salt_1_PA, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,8)
axs[0,1].set_ylim(2500, 0)
#axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,1].tick_params(labelsize=11)
axs[0,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,1].set_title('b) Salinity', fontsize=15)

CS3 = axs[1,0].contourf(lat, depth,  (drho_dy2_PA - drho_dy1_PA)/1e-7, levels = np.linspace(-0.0000005/1e-7, 0.0000005/1e-7, 21), extend = 'both', cmap = 'PuOr_r')
contour1 = axs[1,0].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx0, lat_idx]], colors = 'black')
contour2 = axs[1,0].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx0, lat_idx]], linestyles = '--', colors='black')
contour3 = axs[1,0].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx1, lat_idx]], colors = 'black')
contour4 = axs[1,0].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx1, lat_idx]], linestyles = '--', colors='black')
contour5 = axs[1,0].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx2, lat_idx]], colors = 'black')
contour6 = axs[1,0].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx2, lat_idx]], linestyles = '--', colors='black')

cbar	= colorbar(CS3, ticks = np.arange(-0.0000005/1e-7, 0.00000051/1e-7, 0.0000005/1e-7))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta \rho / \Delta y$ difference [x 10$^-7$ kg/m$^4$]', fontsize = 12)

axs[1,0].set_title('c) Meridional density gradient', fontsize=15)
axs[1,0].set_xlim(-75,8)
axs[1,0].set_ylim(2500, 0)
#axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,0].tick_params(labelsize=11)
axs[1,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,0].set_ylabel('Depth [m]', fontsize=14)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

CS2 = axs[1,1].contourf(lat, depth, UVEL_2_PA - UVEL_1_PA, levels= np.linspace(-0.05, 0.05, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.05, 0.051, 0.05))
cbar.set_label(r'UVEL difference [m/s]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[1,1].set_title('d) Zonal velocity', fontsize=15)
axs[1,1].set_xlim(-75, 8)
axs[1,1].set_ylim(2500, 0)
#axs[1,1].contour(lat, depth, UVEL_1_PA, levels= [0.05], colors = 'black')
axs[1, 1].axvline(x=-65, linestyle='--', color='black')
axs[1, 1].axvline(x=-45, linestyle='--', color='black')
axs[1,1].tick_params(labelsize=11)
axs[1,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
#axs[1,1].set_ylabel('Depth [m]', fontsize=13)
#axs[1,1].tick_params(axis='x', labelsize=11)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

plt.suptitle('Pacific sector', fontsize=19)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_S3_SOM_OS.pdf')
plt.show()

#%%

fig, axs = plt.subplots(3, 2, figsize=(14, 10))

#Temperature
CS = axs[0,0].contourf(lat, depth, temp_2_PA - temp_1_PA, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')

axs[0,0].set_xlim(-75,8)
axs[0,0].set_ylim(2500, 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,0].tick_params(labelsize=12)
#axs[0,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[0,0].set_ylabel('Depth [m]', fontsize=15)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,0].set_title('a) Temperature', fontsize=17)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(temp_1_PA_time[:, depth_i, lat_i], temp_2_PA_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,0].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=15, facecolors="none")

#Salinity
CS2 = axs[0,1].contourf(lat, depth, salt_2_PA - salt_1_PA, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,8)
axs[0,1].set_ylim(2500, 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[0,1].tick_params(labelsize=12)
#axs[0,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 12)
cbar.ax.tick_params(labelsize=12)
axs[0,1].set_title('b) Salinity', fontsize=17)

# non-significance markers for SF (sig < 0.95)
for lat_i in range(0, len(lat), 30):
    for depth_i in range(0, len(depth), 3):
        sig = Welch(salt_1_PA_time[:, depth_i, lat_i], salt_2_PA_time[:, depth_i, lat_i])
        if sig < 0.95:
            axs[0,1].scatter(lat[lat_i], depth[depth_i], marker="o", edgecolor="k", s=15, facecolors="none")


#Density
CS3 = axs[1,0].contourf(lat, depth,  dens_2_PA - dens_1_PA, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

contour1 = axs[1, 0].contour(lat, depth, dens_1_PA, levels=[1027.0], colors='black')
contour2 = axs[1, 0].contour(lat, depth, dens_2_PA, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1, 0].contour(lat, depth, dens_1_PA, levels=[1027.5], colors='black')
contour4 = axs[1,  0].contour(lat, depth, dens_2_PA, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  0].contour(lat, depth, dens_1_PA, levels=[1027.7], colors='black')
contour6 = axs[1,  0].contour(lat, depth, dens_2_PA, levels=[1027.7], linestyles='--', colors='black')

axs[1, 0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.5, 0.51, 0.2))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'Density difference [kg/m$^3$]', fontsize = 12)

axs[1,0].set_title('c) Potential density', fontsize=17)
axs[1,0].set_xlim(-75,8)
axs[1,0].set_ylim(2500, 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,0].tick_params(labelsize=12)
#axs[1,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,0].set_ylabel('Depth [m]', fontsize=15)
#axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=14)

#Meridional density gradient
CS3 = axs[1,1].contourf(lat, depth,  (drho_dy2_PA - drho_dy1_PA)/1e-7, levels = np.linspace(-0.0000005/1e-7, 0.0000005/1e-7, 21), extend = 'both', cmap = 'PuOr_r')

contour1 = axs[1, 1].contour(lat, depth, dens_1_PA, levels=[1027.0], colors='black')
contour2 = axs[1, 1].contour(lat, depth, dens_2_PA, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[1,  1].contour(lat, depth, dens_1_PA, levels=[1027.5], colors='black')
contour4 = axs[1,  1].contour(lat, depth, dens_2_PA, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[1,  1].contour(lat, depth, dens_1_PA, levels=[1027.7], colors='black')
contour6 = axs[1,  1].contour(lat, depth, dens_2_PA, levels=[1027.7], linestyles='--', colors='black')

axs[1, 1].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[1, 1].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005/1e-7, 0.00000051/1e-7, 0.0000005/1e-7))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta \rho / \Delta y$ difference [x 10$^-7$ kg/m$^4$]', fontsize = 12)

axs[1,1].set_title('d) Meridional density gradient', fontsize=17)
axs[1,1].set_xlim(-75,8)
axs[1,1].set_ylim(2500, 0)
axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[1,1].tick_params(labelsize=12)
#axs[1,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[1,1].set_ylabel('Depth [m]', fontsize=15)
#axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=14)


#Stratification frequency
CS3 = axs[2,0].contourf(lat, depth,  N2_lat_2_PA - N2_lat_1_PA, levels= np.linspace(-5e-6, 5e-6, 21), extend = 'both', cmap = 'RdBu_r')

contour1 = axs[2, 0].contour(lat, depth, dens_1_PA, levels=[1027.0], colors='black')
contour2 = axs[2, 0].contour(lat, depth, dens_2_PA, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[2, 0].contour(lat, depth, dens_1_PA, levels=[1027.5], colors='black')
contour4 = axs[2,  0].contour(lat, depth, dens_2_PA, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[2,  0].contour(lat, depth, dens_1_PA, levels=[1027.7], colors='black')
contour6 = axs[2,  0].contour(lat, depth, dens_2_PA, levels=[1027.7], linestyles='--', colors='black')

axs[2, 0].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 0].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-5e-6, 5e-6 + 1e-7, 1e-6))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 12)

axs[2,0].set_title('e) N$^2$', fontsize=17)
axs[2,0].set_xlim(-75,8)
axs[2,0].set_ylim(2500, 0)
axs[2,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[2,0].tick_params(labelsize=12)
#axs[2,0].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[2,0].set_ylabel('Depth [m]', fontsize=15)
axs[2,0].set_xlabel('Latitude [$^\circ$S]', fontsize=15)

#Stratification frequency
CS3 = axs[2,1].contourf(lat, depth,  dN2_S_PA, levels= np.linspace(-5e-6, 5e-6, 21), extend = 'both', cmap = 'RdBu_r')

contour1 = axs[2, 1].contour(lat, depth, dens_1_PA, levels=[1027.0], colors='black')
contour2 = axs[2, 1].contour(lat, depth, dens_2_PA, levels=[1027.0], linestyles='--', colors='black')
contour3 = axs[2, 1].contour(lat, depth, dens_1_PA, levels=[1027.5], colors='black')
contour4 = axs[2,  1].contour(lat, depth, dens_2_PA, levels=[1027.5], linestyles='--', colors='black')
contour5 = axs[2,  1].contour(lat, depth, dens_1_PA, levels=[1027.7], colors='black')
contour6 = axs[2,  1].contour(lat, depth, dens_2_PA, levels=[1027.7], linestyles='--', colors='black')

axs[2, 1].clabel(contour1, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour2, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour3, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour4, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour5, inline=True, fontsize=10, fmt='%1.1f')
axs[2, 1].clabel(contour6, inline=True, fontsize=10, fmt='%1.1f')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-5e-6, 5e-6 + 1e-7, 1e-6))
cbar.ax.tick_params(labelsize=12)
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 12)

axs[2,1].set_title('f) Salinity-driven N$^2$', fontsize=17)
axs[2,1].set_xlim(-75,8)
axs[2,1].set_ylim(2500, 0)
axs[2,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10', '0', '-10'])
axs[2,1].tick_params(labelsize=12)
#axs[2,1].set_xticklabels(['80', '60',  '40',  '20', '0'])
axs[2,1].set_ylabel('Depth [m]', fontsize=15)
axs[2,1].set_xlabel('Latitude [$^\circ$S]', fontsize=15)

plt.suptitle('Pacific sector', fontsize=19)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_S3_SOM_OS_revised.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2_AT - temp_1_AT, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_2_AT - salt_1_AT, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy2_AT - drho_dy1_AT, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')
axs[1,0].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, dens_2_AT - dens_1_AT, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Density', fontsize=14)
axs[1,1].set_xlim(-75, 10)
axs[1,1].set_ylim(depth[-1], 0)
axs[1,1].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx0, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx1, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx2, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx2, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_AT, levels = [dens_1_AT[depth_idx3, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_AT, levels = [dens_2_AT[depth_idx3, lat_idx]], linestyles = '--', colors='black')
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Atlantic sector', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'TEMP_SALT_DENS_Atlantic_OS_HR_POP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2_PA - temp_1_PA, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_2_PA - salt_1_PA, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy2_PA - drho_dy1_PA, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')
axs[1,0].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx0, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx1, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx2, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx2, lat_idx]], linestyles = '--', colors='black')
#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, dens_2_PA - dens_1_PA, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Density', fontsize=14)
axs[1,1].set_xlim(-75, 10)
axs[1,1].set_ylim(depth[-1], 0)
axs[1,1].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx0, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx1, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx2, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx2, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_PA, levels = [dens_1_PA[depth_idx3, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_PA, levels = [dens_2_PA[depth_idx3, lat_idx]], linestyles = '--', colors='black')
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Pacific sector (150$^\circ$E - 295$^\circ$E)', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'TEMP_SALT_DENS_Pacific_OS_HR_POP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2_IN - temp_1_IN, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_2_IN - salt_1_IN, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy2_IN - drho_dy1_IN, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx0, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx1, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx2, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx2, lat_idx]], linestyles = '--', colors='black')
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, dens_2_IN - dens_1_IN, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Density', fontsize=14)
axs[1,1].set_xlim(-75, 10)
axs[1,1].set_ylim(depth[-1], 0)
axs[1,1].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx0, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx1, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx2, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx2, lat_idx]], linestyles = '--', colors='black')
axs[1,1].contour(lat, depth, dens_1_IN, levels = [dens_1_IN[depth_idx3, lat_idx]], colors = 'black')
axs[1,1].contour(lat, depth, dens_2_IN, levels = [dens_2_IN[depth_idx3, lat_idx]], linestyles = '--', colors='black')
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Indian sector (25$^\circ$E - 150$^\circ$E)', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'TEMP_SALT_DENS_Indian_OS_HR_POP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, dR_T_AT, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Temperature-driven density change [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature-driven density change', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, dR_S_AT, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity-driven density change', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, dR_lin_AT, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3)
#cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Linearised density anomaly', fontsize=14)
axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, recon_error_AT, levels= np.linspace(-0.02, 0.02, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.02, 0.021, 0.005))
#cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Reconstruction error', fontsize=14)
axs[1,1].set_xlim(-75, 10)
axs[1,1].set_ylim(depth[-1], 0)
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Atlantic sector', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'DENS_components_Atlantic_OS_HR_POP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, dR_T_PA, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Temperature-driven density change [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature-driven density change', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, dR_S_PA, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity-driven density change', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, dR_lin_PA, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3)
#cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Linearised density anomaly', fontsize=14)
axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, recon_error_PA, levels= np.linspace(-0.02, 0.02, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.02, 0.021, 0.005))
#cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Reconstruction error', fontsize=14)
axs[1,1].set_xlim(-75, 10)
axs[1,1].set_ylim(depth[-1], 0)
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Pacific sector', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'DENS_components_Pacific_OS_HR_POP.pdf')
plt.show()

#%%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, dR_T_IN, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Temperature-driven density change [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature-driven density change', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, dR_S_IN, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity-driven density change', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, dR_lin_IN, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3)
#cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Linearised density anomaly', fontsize=14)
axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, recon_error_IN, levels= np.linspace(-0.02, 0.02, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.02, 0.021, 0.005))
#cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Reconstruction error', fontsize=14)
axs[1,1].set_xlim(-75, 10)
axs[1,1].set_ylim(depth[-1], 0)
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Indian sector', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'DENS_components_Indian_OS_HR_POP.pdf')
plt.show()



































# %%
