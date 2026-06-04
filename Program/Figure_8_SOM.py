#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 22 15:28:04 2025

@author: 6008399

Figure 8 SOM OS

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

time_mld        = fh.variables['time'][:]       #time [model years]
MLD_SO_PA       = fh.variables['MLD'][:, 200:600,0:500]  
lat         = fh.variables['lat'][:]    #latitude
lon         = fh.variables['lon'][:]    #latitude
lon_PA          = lon[0:500]
lat_PA          = lat[200:600]

fh.close()

fh = netcdf.Dataset(directory_zenodo + 'MLD_year_1-600_NZ_160_190E_70S_62S.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_NZ = fh.variables['MLD_max'][:] 
MLD = fh.variables['MLD'][:]

fh.close()

fh = netcdf.Dataset(directory_zenodo + 'MLD_year_1-600_AU_70S_50S.nc','r')

time_mld               = fh.variables['time'][:]       #time [model years]
MLD_AU             = fh.variables['MLD'][:]        #Transport [Sv]
MLD_max_AU         = fh.variables['MLD_max'][:]
lon_AU             = fh.variables['lon'][:]        #Longitude
lat_AU             = fh.variables['lat'][:]        #Latitude

fh.close()

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld        = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_WGKP    = fh.variables['MLD_max'][:] 

fh.close()

#%%

print(lat_PA[0], lat_PA[-1])
print(lon_PA[0], lon_PA[-1])

#%%

N2_time_WGKP    = ma.masked_all((len(time_all), len(depth_PD)))
N2_time_NZ      = ma.masked_all((len(time_all), len(depth_PD)))
N2_time_AU      = ma.masked_all((len(time_all), len(depth_PD)))
N2_time_PA      = ma.masked_all((len(time_all), len(depth_PD)))

for time_i in range(len(time_all)):
    n0, N2_time_PA[time_i, :] = compute_N2_from_profile(PD_PA_all[time_i,:], depth_PD)
    n0, N2_time_AU[time_i, :] = compute_N2_from_profile(PD_AU_all[time_i,:], depth_PD)
    n0, N2_time_NZ[time_i, :] = compute_N2_from_profile(PD_NZ_all[time_i,:], depth_PD)
    n0, N2_time_WGKP[time_i, :] = compute_N2_from_profile(PD_WGKP_all[time_i,:], depth_PD)    
    
#%%

window = 20 

fig, ax = plt.subplots(figsize=(8, 5)) 
CS = ax.contourf(
    time_all, depth_PD, (N2_time_AU - np.mean(N2_time_AU[1:100], axis=0)).transpose(),
    levels=np.linspace(-0.0000006, 0.0000006, 21), extend='both', cmap='RdBu_r'
)

ax.plot(
    time_mld[window//2:-window//2+1], Moving_average(MLD_max_AU, window), linewidth=2,
    color='black', label='maximum MLD'
)

ax.set_title('c) AU convective region', fontsize=14)
ax.set_ylim(2000, 1)  # Invert y-axis
ax.set_xlim(1, 600)
ax.set_ylabel('Depth [m]', fontsize=12)
ax.set_xlabel('Time [model years]', fontsize=12)
ax.legend(loc=4, fontsize=11)
cbar = fig.colorbar(CS, ax=ax)#, ticks=[-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])  # Set specific ticks
cbar.set_label('Area-averaged $N^2$ anomaly [s$^{-2}$]', fontsize=11)
plt.tight_layout()
plt.savefig(directory_figures + 'AU_N2_MLD.pdf')
    
#%%

window = 20 

fig, ax = plt.subplots(figsize=(8, 5)) 
CS = ax.contourf(
    time_all, depth_PD, (N2_time_WGKP - np.mean(N2_time_WGKP[1:100], axis=0)).transpose(),
    levels=np.linspace(-0.0000006, 0.0000006, 21), extend='both', cmap='RdBu_r'
)

ax.plot(
    time_mld[window//2:-window//2+1], Moving_average(MLD_max_WGKP, window), linewidth=2,
    color='black', label='maximum MLD'
)

ax.set_title('b) WGKP convective region', fontsize=14)
ax.set_ylim(2000, 1)  # Invert y-axis
ax.set_xlim(1, 600)
ax.set_ylabel('Depth [m]', fontsize=12)
ax.set_xlabel('Time [model years]', fontsize=12)
ax.legend(loc=4, fontsize=11)
cbar = fig.colorbar(CS, ax=ax)#, ticks=[-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])  # Set specific ticks
cbar.set_label('Area-averaged $N^2$ anomaly [s$^{-2}$]', fontsize=11)
plt.tight_layout()
plt.savefig(directory_figures + 'WGKP_N2_MLD.pdf')
    
#%%

window = 20

fig, ax = plt.subplots(figsize=(8, 5)) 
CS = ax.contourf(
    time_all, depth_PD, (N2_time_NZ - np.mean(N2_time_NZ[1:100], axis=0)).transpose(),
    levels=np.linspace(-0.0000006, 0.0000006, 21), extend='both', cmap='RdBu_r'
)

ax.plot(
    time_mld[window//2:-window//2+1], Moving_average(MLD_max_NZ, window), linewidth=2,
    color='black', label='maximum MLD'
)

ax.set_title('a) NZ convective region', fontsize=14)
ax.set_ylim(2000, 1)  # Invert y-axis
ax.set_xlim(1, 600)
ax.set_ylabel('Depth [m]', fontsize=12)
ax.set_xlabel('Time [model years]', fontsize=12)
ax.legend(loc=4, fontsize=11)
cbar = fig.colorbar(CS, ax=ax)#, ticks=[-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])  # Set specific ticks
cbar.set_label('Area-averaged $N^2$ anomaly [s$^{-2}$]', fontsize=11)
plt.tight_layout()
plt.savefig(directory_figures + 'NZ_N2_MLD.pdf')
plt.show()
    
    
#%%

window = 20 

fig, ax = plt.subplots(figsize=(8, 5)) 
CS = ax.contourf(
    time_all, depth_PD, (N2_time_PA - np.mean(N2_time_PA[1:100], axis=0)).transpose(),
    levels=np.linspace(-0.0000006, 0.0000006, 21), extend='both', cmap='RdBu_r'
)

ax.plot(
    time_mld[window//2:-window//2+1], Moving_average(np.max(MLD_SO_PA, axis=(1, 2)), window), linewidth=2,
    color='black', label='maximum MLD'
)

ax.set_title('d) PA convective region', fontsize=14)
ax.set_ylim(2000, 1)  # Invert y-axis
ax.set_xlim(1, 600)
ax.set_ylabel('Depth [m]', fontsize=12)
ax.set_xlabel('Time [model years]', fontsize=12)
ax.legend(loc=4, fontsize=11)
cbar = fig.colorbar(CS, ax=ax)#, ticks=[-0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06])  # Set specific ticks
cbar.set_label('Area-averaged $N^2$ anomaly [s$^{-2}$]', fontsize=11)
plt.tight_layout()
plt.savefig(directory_figures + 'PA_N2_MLD.pdf')

#%%

PD_profiles = {
    "NZ": PD_NZ_all,
    "WGKP": PD_WGKP_all,
    "AU": PD_AU_all,
    "PA": PD_PA_all,
}

MLD_max = {
    "NZ": MLD_max_NZ,
    "WGKP": MLD_max_WGKP,
    "AU": MLD_max_AU,
    "PA": np.max(MLD_SO_PA, axis=(1, 2)),
}

N2_time_all = {}
for region, profile in PD_profiles.items():
    N2_time = ma.masked_all((len(time_all), len(depth_PD)))
    for time_i in range(len(time_all)):
        _, N2_time[time_i, :] = compute_N2_from_profile(profile[time_i, :], depth_PD)
    N2_time_all[region] = N2_time

top_levels = np.linspace(-0.00001/1e-6, 0.00001/1e-6, 21) 
bot_levels = np.linspace(-0.00002/1e-6, 0.00002/1e-6, 21) 

fig, axs = plt.subplots(2, 2, figsize=(10, 6))

titles = { "NZ": "a) NZ convective region", "WGKP": "b) WGKP convective region", "AU": "c) AU convective region", "PA": "d) PA convective region", }

for i, (ax, (region, N2_time)) in enumerate(zip(axs.flat, N2_time_all.items())): 
    data = ((N2_time - np.mean(N2_time[1:100], axis=0))/1e-7).T 
    levels = top_levels if i < 2 else bot_levels
    CS = ax.contourf(
    time_all, depth_PD, data,
    levels=levels, extend='both', cmap='RdBu_r')

    ax.plot(
    time_all[window // 2:-window // 2 + 1],
    np.convolve(MLD_max[region], np.ones(window)/window, mode='valid'),
    linewidth=2, color='black', label='maximum MLD')
    ax.set_title(titles[region], fontsize=14)
    ax.set_ylim(2000, -1)
    ax.set_xlim(-1, 600)
    ax.legend(loc=4, fontsize=11)

    # Colorbar per subplot using the same levels range as the artist
    cbar = fig.colorbar(CS, ax=ax)
    cbar.set_label('N$^2$ anomaly [x 10$^{-7}$ s$^{-2}$]', fontsize=11)

axs[0,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Time [model year]', fontsize=12)
axs[1,1].set_xlabel('Time [model year]', fontsize=12)
#for ax in axs.flat: 
#    ax.label_outer() # hides x on top row and y on right column

#fig.supxlabel('Time [model years]', fontsize=12) 
#fig.supylabel('Depth [m]', fontsize=12)

plt.tight_layout() 
plt.savefig(directory_figures +'Figure_8_SOM_OS.pdf')
plt.show()


# %%
