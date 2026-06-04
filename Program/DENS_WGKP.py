#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 16 09:29:52 2025

@author: 6008399

Density in WGKP region

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
directory_data = r'/Users/6008399/Documents/PhD/HR_POP/Zenodo/Data_Final/'
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
    
    if N2.any() < 0:
        print('N2 is unstable!')

    return n0, N2

# def compute_N2_from_profile(density, depth, rho0=1027.0, g=9.81):
#     density = np.asarray(density)
#     depth = np.asarray(depth)

#     drho_dz = np.gradient(density, depth)
#     N2 = (g / rho0) * drho_dz   # because depth increases downward

#     if np.any(N2 < 0):
#         print("N2 is unstable!")

#     return drho_dz, N2

def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#%% Mixed layer depth WGKP

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld        = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_WGKP    = fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Volume_90_50S_35W_80E_WGKP.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat_vol         = fh.variables['lat'][:]    #latitude
lon_vol         = fh.variables['lon'][:]    #latitude
volume      = fh.variables['volume'][:]   #salinity

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/PD_gsw_year_1-600_area_averaged_WGKP.nc','r')

depth_PD       = fh.variables['depth'][:]  #depth
time_all = fh.variables['time'][:] #Time SOM cycle 1
PD_WGKP_all = fh.variables['PD'][:]
PD_WGKP_1 = fh.variables['PD'][1:100]
PD_WGKP_2 = fh.variables['PD'][500:600]

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_zonal_averaged_-35E-80E_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1      = fh.variables['TEMP'][:]   #temperature
salt_1      = fh.variables['SALT'][:]   #salinity
dens_1      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_300-400_zonal_averaged_-35E-80E_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_300      = fh.variables['TEMP'][:]   #temperature
salt_300      = fh.variables['SALT'][:]   #salinity
dens_300      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_400-500_zonal_averaged_-35E-80E_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_400      = fh.variables['TEMP'][:]   #temperature
salt_400      = fh.variables['SALT'][:]   #salinity
dens_400      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_-35E-80E_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_500      = fh.variables['TEMP'][:]   #temperature
salt_500      = fh.variables['SALT'][:]   #salinity
dens_500      = fh.variables['PD'][:]     #potential density

fh.close()

#%% Meridional density gradient

bin_size = 20
drho_dy1 = ma.masked_all((len(depth), len(lat)))
drho_dy300 = ma.masked_all((len(depth), len(lat)))
drho_dy400 = ma.masked_all((len(depth), len(lat)))
drho_dy500 = ma.masked_all((len(depth), len(lat)))

deg2m = 111e3   # degrees latitude → meters

for depth_i in range(len(depth)):
    for lat_i in range(bin_size, len(lat) - bin_size):

        lat_1 = lat_i - bin_size
        lat_2 = lat_i + bin_size
        
        print(lat_1)
        print(lat_2)
        
        #sys.exit()

        # actual latitude difference (in meters)
        del_lat_m = np.abs(lat[lat_2] - lat[lat_1])  * deg2m

        # density difference
        del_dens1 = dens_1[depth_i, lat_2] - dens_1[depth_i, lat_1]
        del_dens300 = dens_300[depth_i, lat_2] - dens_300[depth_i, lat_1]
        del_dens400 = dens_400[depth_i, lat_2] - dens_400[depth_i, lat_1]
        del_dens500 = dens_500[depth_i, lat_2] - dens_500[depth_i, lat_1]

        drho_dy1[depth_i, lat_i] = del_dens1 / del_lat_m
        drho_dy300[depth_i, lat_i] = del_dens300 / del_lat_m
        drho_dy400[depth_i, lat_i] = del_dens400 / del_lat_m
        drho_dy500[depth_i, lat_i] = del_dens500 / del_lat_m

#%% WGKP densities

plt.figure()
plt.plot(np.mean(PD_WGKP_1, axis=0), depth)
plt.plot(np.mean(dens_1, axis=1), depth, label='before')

plt.plot(np.mean(PD_WGKP_2, axis=0), depth)
plt.plot(np.mean(dens_500, axis=1), depth, label='after')
plt.legend()
plt.ylim(2000, 0)

#%% WGKP densities

volume_lat = np.sum(volume, axis=2)

lat_idx = np.where(lat == lat_vol[0])[0][0] 
lat_idx_end = np.where(lat == lat_vol[-1])[0][0]

print(lat[lat_idx])
print(lat[lat_idx_end])

PD_weighted_1 = np.sum(dens_1[:, lat_idx:lat_idx_end+1] * volume_lat, axis=1) / np.sum(volume_lat, axis=1)
PD_weighted_2 = np.sum(dens_500[:, lat_idx:lat_idx_end+1] * volume_lat, axis=1) / np.sum(volume_lat, axis=1)

plt.figure()
plt.plot(PD_weighted_1, depth, label='before')
plt.plot(PD_weighted_2, depth, label='after')
plt.legend()
plt.ylim(2000, 0)

#%%Test

n0 = np.zeros(len(depth))
for depth_i in range(len(depth)):
	if depth_i == 0:
		n0[depth_i] = (PD_weighted_1[depth_i] - PD_weighted_1[depth_i+1])/(depth[1] - depth[0])
	elif depth_i == len(depth) - 1:
		n0[depth_i] = (PD_weighted_1[depth_i - 1] - PD_weighted_1[depth_i])/(depth[-1] - depth[-2])
	else:
		n0[depth_i] = (PD_weighted_1[depth_i-1] - PD_weighted_1[depth_i+1])/(np.abs(depth[depth_i-1] - depth[depth_i+1]))

#%%

rho0 = 1027
g = 9.81

N2 = - g/rho0 * n0

#Something weird happens in the lower layers (PD profiles are also reaally weird here. Upper 3000m is correct though)
plt.figure()
plt.plot(N2, depth)

#%%

n0_som_1_WGKP, N2_som_1_WGKP_test = compute_N2_from_profile(PD_weighted_1, depth)
n0_som_2_WGKP, N2_som_2_WGKP = compute_N2_from_profile(PD_weighted_2, depth)

n0_som_1_WGKP, N2_som_1_WGKP = compute_N2_from_profile(np.mean(PD_WGKP_1, axis = 0), depth)
n0_som_2_WGKP, N2_som_2_WGKP = compute_N2_from_profile(np.mean(PD_WGKP_2, axis=0), depth)

plt.figure()
plt.plot(N2_som_1_WGKP, depth)
plt.plot(N2_som_1_WGKP_test, depth)
plt.ylim(200, 0)

#%%
n0_lat_1 = ma.masked_all((len(depth), len(lat)))
n0_lat_300 = ma.masked_all((len(depth), len(lat)))
n0_lat_400 = ma.masked_all((len(depth), len(lat)))
n0_lat_500 = ma.masked_all((len(depth), len(lat)))

N2_lat_1 = ma.masked_all((len(depth), len(lat)))
N2_lat_300 = ma.masked_all((len(depth), len(lat)))
N2_lat_400 = ma.masked_all((len(depth), len(lat)))
N2_lat_500 = ma.masked_all((len(depth), len(lat)))

for lat_i in range(len(lat)):
    n0_lat_1[:, lat_i], N2_lat_1[:,lat_i] = compute_N2_from_profile(dens_1[:,lat_i], depth)
    n0_lat_300[:, lat_i], N2_lat_300[:,lat_i] = compute_N2_from_profile(dens_300[:,lat_i], depth)
    n0_lat_400[:, lat_i], N2_lat_400[:,lat_i] = compute_N2_from_profile(dens_400[:,lat_i], depth)
    n0_lat_500[:, lat_i], N2_lat_500[:,lat_i] = compute_N2_from_profile(dens_500[:,lat_i], depth)

plt.figure()
plt.plot(N2_som_1_WGKP, depth)
plt.plot(N2, depth)
plt.plot(np.nanmean(N2_lat_1, axis=1), depth, '--', label='N2 from lat profiles') #Not area-averged, and also mean over 80S - 10N!!
#plt.plot(np.nanmean(N2_lat_300, axis=1), depth)
#plt.plot(np.mean(N2_lat_400, axis=1), depth)
#plt.plot(np.mean(N2_lat_500, axis=1), depth)
plt.ylim(2000, 0)

N2_time = ma.masked_all((len(time_all), len(depth)))
for time_i in range(len(time_all)):
    n0, N2_time[time_i, :] = compute_N2_from_profile(PD_WGKP_all[time_i,:], depth)
    
plt.figure()
plt.contourf(lat, depth, N2_lat_1)
    
#%%

window = 20 

fig, ax = plt.subplots(figsize=(8, 5)) 
CS = ax.contourf(
    time_all, depth_PD, (N2_time - np.mean(N2_time[1:100], axis=0)).transpose(),
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

#%% Isopycnal slope

s1 = - (g/rho0) * drho_dy1/N2_lat_1
s2 = - (g/rho0) * drho_dy500/N2_lat_500

#s1[-10:, :] = np.ma.masked
#s2[-10:, :] = np.ma.masked

#%%

lat_idx = -1 #latitude of 0.05 degN (equator)
depth_idx0 = 18 #depth of 465m
depth_idx1 = 21 #depth of 918m
depth_idx2 = 25 #depth of 1380m
depth_idx3 = 35 #depth of 2125m

plt.figure(figsize=(12,6))
plt.title('Isopycnal slope')
plt.contourf(lat, depth, s2 - s1, levels = np.arange(-0.01, 0.01, 0.001), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [dens_500[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [dens_500[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [dens_500[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [dens_500[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1], 0)
#plt.xlim(-80, -50)

#%%

window = 5

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS3 = ax1.contourf(lat, depth, drho_dy500 - drho_dy1, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
ax1.set_title('a) Merdidional density gradient', fontsize=14)
ax1.set_xlim(-80, -50)
ax1.set_ylim(depth[-1], 0)
#ax1.set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = ax2.contourf(lat, depth, dens_500 - dens_1, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
ax2.set_title('b) Density difference', fontsize=14)
ax2.set_xlim(-80, -50)
ax2.set_ylim(depth[-1], 0)
#ax2.contour(lat, depth, dens_1, levels = [dens_1[depth_idx0, lat_idx]], colors = 'black')
#ax2.contour(lat, depth, dens_2, levels = [dens_2[depth_idx0, lat_idx]], linestyles = '--', colors='black')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS = ax3.contourf(time_all, depth_PD, (PD_WGKP_all - np.mean(PD_WGKP_all[1:100], axis=0)).transpose(), levels=np.linspace(-0.06, .06, 21), extend='both', cmap = 'RdBu_r')
cbar = fig.colorbar(CS, ax=ax3, ticks=[-0.06, -0.04, -0.02,  0,  0.02,  0.04, 0.06])  # Set specific ticks
cbar.set_label('Potential density anomaly [kg/m$^3$]', fontsize = 11)
ax3.plot(time_mld[window//2:-window//2+1], Moving_average(MLD_max_WGKP, window), color='black', label='maximum MLD')
ax3.set_title('c) MLD maximum and PD anomaly', fontsize=14)
#ax1.plot(time_mld, np.mean(MLD, axis=(1,2)))
ax3.set_ylim(2000, 1)
ax3.set_xlim(1,600)
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Time [model years]', fontsize=12)
ax3.legend(loc=4, fontsize=11)

#ax4.plot(N2_som_2_WGKP - N2_som_1_WGKP, depth, color='blue', label='SOM2 vs SOM1')
ax4.plot(N2_som_2_WGKP - N2_som_1_WGKP, depth, color='red', label='(500-600) vs (1-100)')
ax4.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax4.set_title('d) Area-averaged N$^2$ difference', fontsize=14)
ax4.set_xlabel('N$^2$ difference [s$^{-1}$]', fontsize=12)
#ax2.set_ylabel('Depth [m]')
ax4.legend(loc=4, fontsize=11)
ax4.set_ylim(5000, 1)
ax4.set_ylim(2000, 1)
ax4.set_xlim(-4.e-6, 10e-6)

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$E, 80 - 50$^\circ$S)', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Dens_MLD_N2_WGKP.pdf')



fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lat, depth, dens_1, levels= np.linspace(1025, 1027.7, 21), extend = 'both', cmap = 'Spectral_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(1025, 1027.8, 0.5))
cbar.set_label(r'$\rho$ [kg/m$^4$]', fontsize = 11)
ax1.set_title('a) Density (1 - 100)', fontsize=14)
ax1.set_xlim(-80,-30)
ax1.set_ylim(2000, 0)
cs = ax1.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.1], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax1.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax1.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax2.contourf(lat, depth, dens_300 - dens_1, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\Delta \rho$ [kg/m$^4$]', fontsize = 11)
ax2.set_title('b) Density difference (300-400) vs (1 - 100)', fontsize=14)
ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
cs = ax2.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.0], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.7], linestyles = '--', colors='black')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lat, depth, dens_400 - dens_1, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\Delta \rho$ [kg/m$^4$]', fontsize = 11)
ax3.set_title('c) Density difference (400-500) vs (1-100)', fontsize=14)
ax3.set_xlim(-80,-30)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.0], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.5], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.7], linestyles = '--', colors='black')
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax3.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax3.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax4.contourf(lat, depth, dens_500 - dens_1, levels= np.linspace(-0.2, 0.2, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.2, 0.21, 0.05))
cbar.set_label(r'$\Delta \rho$ [kg/m$^4$]', fontsize = 11)
ax4.set_title('d) Density difference (500-600) vs (1-100)', fontsize=14)
ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
#ax4.contour(lat, depth, dens_1, levels = [dens_1[depth_idx0, lat_idx]], colors = 'black')
#ax4.contour(lat, depth, dens_2, levels = [dens_2[depth_idx0, lat_idx]], linestyles = '--', colors='black')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax4.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax4.vlines(x=-50, ymin=2000, ymax=0, color='grey')

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Density_evolving_WGKP.pdf')

# %%

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lat, depth, drho_dy1, levels= np.linspace(-2.6e-6, 2.6e-6, 21), extend = 'both', cmap = 'seismic_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-2.4e-6, 2.8e-6, 1e-6))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
ax1.set_title('a) Meridional density gradient (1 - 100)', fontsize=14)
ax1.set_xlim(-80,-30)
ax1.set_ylim(2000, 0)
cs = ax1.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.1], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax1.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax1.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax2.contourf(lat, depth, drho_dy300 - drho_dy1, levels= np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.0000005, 0.00000051, 0.0000001))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
ax2.set_title('b) Meridional density gradient difference (300-400) vs (1 - 100)', fontsize=14)
ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
cs = ax2.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.0], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.7], linestyles = '--', colors='black')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lat, depth, drho_dy400 - drho_dy1, levels= np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.0000005, 0.00000051, 0.0000001))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
ax3.set_title('c) Meridional density gradient difference (400-500) vs (1-100)', fontsize=14)
ax3.set_xlim(-80,-30)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.0], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.5], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.7], linestyles = '--', colors='black')
#ax1.set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax3.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax3.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax4.contourf(lat, depth, drho_dy500 - drho_dy1, levels= np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.0000005, 0.00000051, 0.0000001))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
ax4.set_title('d) Meridional density gradient difference (500-600) vs (1-100)', fontsize=14)
ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
#ax4.contour(lat, depth, dens_1, levels = [dens_1[depth_idx0, lat_idx]], colors = 'black')
#ax4.contour(lat, depth, dens_2, levels = [dens_2[depth_idx0, lat_idx]], linestyles = '--', colors='black')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax4.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax4.vlines(x=-50, ymin=2000, ymax=0, color='grey')

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Meridional_density_gradients_WGKP.pdf')

# %%
# %%

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lat, depth, salt_1, levels= np.linspace(33, 35, 21), extend = 'both')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(33, 35.1, 0.5))
cbar.set_label(r'Salinity [g/kg]', fontsize = 11)
ax1.set_title('a) Salinity (1 - 100)', fontsize=14)
ax1.set_xlim(-80,-30)
ax1.set_ylim(2000, 0)
cs = ax1.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.1], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax1.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax1.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax2.contourf(lat, depth, salt_300 - salt_1, levels= np.linspace(-1, 1, 21), extend = 'both', cmap = 'BrBG_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.1, 0.2))
cbar.set_label(r'$\Delta S$ [g/kg]', fontsize = 11)
ax2.set_title('b) Salinity difference (300-400) vs (1 - 100)', fontsize=14)
ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
cs = ax2.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.0], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.7], linestyles = '--', colors='black')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lat, depth, salt_400 - salt_1, levels= np.linspace(-1, 1, 21), extend = 'both', cmap = 'BrBG_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.1, 0.2))
cbar.set_label(r'$\Delta S$ [g/kg]', fontsize = 11)
ax3.set_title('c) Salinity difference (400-500) vs (1-100)', fontsize=14)
ax3.set_xlim(-80,-30)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.0], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.5], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.7], linestyles = '--', colors='black')
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax3.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax3.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax4.contourf(lat, depth, salt_500 - salt_1, levels= np.linspace(-1, 1, 21), extend = 'both', cmap = 'BrBG_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.1, 0.2))
cbar.set_label(r'$\Delta S$ [g/kg]', fontsize = 11)
ax4.set_title('d) Salinity difference (500-600) vs (1-100)', fontsize=14)
ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
cs = ax4.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax4.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax4.vlines(x=-50, ymin=2000, ymax=0, color='grey')

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Salinities_WGKP.pdf')

# %%

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lat, depth, temp_1, levels= np.linspace(0, 30, 21), extend = 'both', cmap = 'Spectral_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(0, 30.1, 2))
cbar.set_label(r'Temperature [°C]', fontsize = 11)
ax1.set_title('a) Temperature (1 - 100)', fontsize=14)
ax1.set_xlim(-80,-30)
ax1.set_ylim(2000, 0)
cs = ax1.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.1], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax1.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax1.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax2.contourf(lat, depth, temp_300 - temp_1, levels= np.linspace(-4, 4, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4, 4.1, 1))
cbar.set_label(r'$\Delta T$ [°C]', fontsize = 11)
ax2.set_title('b) Temperature difference (300-400) vs (1 - 100)', fontsize=14)
ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
cs = ax2.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.0], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.7], linestyles = '--', colors='black')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lat, depth, temp_400 - temp_1, levels= np.linspace(-4, 4, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4, 4.1, 1))
cbar.set_label(r'$\Delta T$ [°C]', fontsize = 11)
ax3.set_title('c) Temperature difference (400-500) vs (1-100)', fontsize=14)
ax3.set_xlim(-80,-30)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.0], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.5], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.7], linestyles = '--', colors='black')
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax3.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax3.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax4.contourf(lat, depth, temp_500 - temp_1, levels= np.linspace(-4, 4, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4, 4.1, 1))
cbar.set_label(r'$\Delta T$ [°C]', fontsize = 11)
ax4.set_title('d) Temperature difference (500-600) vs (1-100)', fontsize=14)
ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
cs = ax4.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax4.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax4.vlines(x=-50, ymin=2000, ymax=0, color='grey')

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Temperatures_WGKP.pdf')

# %%
fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lat, depth, N2_lat_1, levels= np.linspace(0, 3.5e-5, 21), extend = 'max', cmap = 'Spectral_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(0, 3.5e-5 + 1e-6, 1e-5))
cbar.set_label(r'N² [s⁻²]', fontsize = 11)
ax1.set_title('a) Buoyancy Frequency (1 - 100)', fontsize=14)
ax1.set_xlim(-80,-30)
ax1.set_ylim(2000, 0)
cs = ax1.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.1], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#ax1.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
cs = ax1.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax1.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax1.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax2.contourf(lat, depth, N2_lat_300 - N2_lat_1, levels= np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-2e-6, 2e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax2.set_title('b) Buoyancy Frequency difference (300-400) vs (1 - 100)', fontsize=14)
ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
cs = ax2.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.0], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.5], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_300, levels = [1027.7], linestyles = '--', colors='black')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lat, depth, N2_lat_400 - N2_lat_1, levels= np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-2e-6, 2e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax3.set_title('c) Buoyancy Frequency difference (400-500) vs (1-100)', fontsize=14)
ax3.set_xlim(-80,-30)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.0], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.5], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_400, levels = [1027.7], linestyles = '--', colors='black')
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax3.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax3.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax4.contourf(lat, depth, N2_lat_500 - N2_lat_1, levels= np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-2e-6, 2e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax4.set_title('d) Buoyancy Frequency difference (500-600) vs (1-100)', fontsize=14)
ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
cs = ax4.contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax4.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax4.vlines(x=-50, ymin=2000, ymax=0, color='grey')

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Buoyancy_frequency_WGKP.pdf')

# %% Temperature and salinity driven components of N2 gradients

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
dT = temp_500 - temp_1
dS = salt_500 - salt_1
dR = dens_500 - dens_1

#Reference state for the derivatives (midpoint state)
T_ref = 0.5 * (temp_1 + temp_500)
S_ref = 0.5 * (salt_1 + salt_500)

drho_dT = RHO_0_dT(T_ref, S_ref)   # shape: (depth, lat)
drho_dS = RHO_0_dS(T_ref, S_ref)

#Decompose density change
dR_T = drho_dT * dT    # temperature-driven density change
dR_S = drho_dS * dS    # salinity-driven density change

#linear reconstructed total density change
dR_lin = dR_T + dR_S

#Reconstruction error
recon_error = dR - dR_lin

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, dR_T, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,0].set_xlim(-75,-30)
axs[0,0].set_ylim(depth[-1], 0)
#axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Temperature-driven density change [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature-driven density change', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, dR_S, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'PuOr_r')
axs[0,1].set_xlim(-75,-30)
axs[0,1].set_ylim(depth[-1], 0)
#axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
#cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity-driven density change', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, dR_lin, levels = np.linspace(-0.5, 0.5, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3)
#cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Linearised density anomaly', fontsize=14)
axs[1,0].set_xlim(-75,-30)
axs[1,0].set_ylim(depth[-1], 0)
#axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS2 = axs[1,1].contourf(lat, depth, recon_error, levels= np.linspace(-0.02, 0.02, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.02, 0.021, 0.005))
#cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Reconstruction error', fontsize=14)
axs[1,1].set_xlim(-75, -30)
axs[1,1].set_ylim(depth[-1], 0)
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'DENS_components_WGKP_OS_HR_POP.pdf')
plt.show()
# %%

g = 9.81
rho0 = 1027.0

def d_dz(var, depth):
    return np.gradient(var, depth, axis=0)

# Actual N² in each state for depth increasing downward!! (so positive N² means stable stratification)
N2_1 = (g/rho0) * d_dz(dens_1, depth)
N2_2 = (g/rho0) * d_dz(dens_500, depth)
dN2 = N2_2 - N2_1

# T/S contributions to ΔN²
dN2_T = (g/rho0) * d_dz(dR_T, depth)
dN2_S = (g/rho0) * d_dz(dR_S, depth)
dN2_lin = dN2_T + dN2_S

# Reconstruction error
recon_error_N2 = dN2 - dN2_lin

# %%

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, dN2_T, levels = np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')
axs[0,0].set_xlim(-75,-30)
axs[0,0].set_ylim(2000, 0)
axs[0,0].contour(lat, depth, dens_1, levels = [1027.1], colors = 'black')
axs[0,0].contour(lat, depth, dens_500, levels = [1027.1], linestyles = '--', colors='black')
axs[0,0].contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
axs[0,0].contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
axs[0,0].contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
axs[0,0].contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
#axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-1e-5, 1.01e-5, 0.5e-5))
#cbar.set_label('Temperature-driven density change [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature-driven N$^2$ change', fontsize=14)
axs[0,0].vlines(x=-80, ymin=2000, ymax=0, color='grey')
axs[0,0].vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = axs[0,1].contourf(lat, depth, dN2_S, levels = np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')
axs[0,1].set_xlim(-75,-30)
axs[0,1].set_ylim(2000, 0)
cs = axs[0,1].contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
axs[0,1].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[0,1].contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
axs[0,1].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[0,1].contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
axs[0,1].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[0,1].contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
axs[0,1].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[0,1].contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
axs[0,1].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[0,1].contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
axs[0,1].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1e-5, 1.01e-5, 0.5e-5))
#cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity-driven N$^2$ change', fontsize=14)
axs[0,1].vlines(x=-80, ymin=2000, ymax=0, color='grey')
axs[0,1].vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS3 = axs[1,0].contourf(lat, depth, dN2_lin, levels = np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS3)
#cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Linearised N$^2$ anomaly', fontsize=14)
axs[1,0].set_xlim(-75,-30)
axs[1,0].set_ylim(2000, 0)
cs = axs[1,0].contour(lat, depth, dens_1, levels = [1027.0], colors = 'black')
axs[1,0].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[1,0].contour(lat, depth, dens_500, levels = [1027.0], linestyles = '--', colors='black')
axs[1,0].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[1,0].contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
axs[1,0].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[1,0].contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
axs[1,0].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[1,0].contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
axs[1,0].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = axs[1,0].contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
axs[1,0].clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
#axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)
axs[1,0].vlines(x=-80, ymin=2000, ymax=0, color='grey')
axs[1,0].vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = axs[1,1].contourf(lat, depth, recon_error_N2, levels= np.linspace(-0.02e-5, 0.02e-5, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-0.02e-5, 0.02e-5, 0.005e-5))
#cbar.set_label(r'$\rho$ difference [kg/m$^3$]', fontsize = 11)
axs[1,1].set_title('d) Reconstruction error', fontsize=14)
axs[1,1].set_xlim(-75, -30)
axs[1,1].set_ylim(2000, 0)
#axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,1].set_ylabel('Depth [m]', fontsize=12)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)
axs[1,1].vlines(x=-80, ymin=2000, ymax=0, color='grey')
axs[1,1].vlines(x=-50, ymin=2000, ymax=0, color='grey')

plt.suptitle('Convection region WGKP (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=16)

plt.tight_layout()
plt.savefig(directory_figures +'N2_components_WGKP_OS_HR_POP.pdf')
plt.show()
# %% Check wheter N2 latitude is consistent with weighte N2 derived from area-averaged PD profile

lat_idx = np.where(lat == lat_vol[0])[0][0] 
lat_idx_end = np.where(lat == lat_vol[-1])[0][0]

plt.figure(figsize=(8, 6))
plt.plot(np.nansum(N2_lat_1[:, lat_idx:lat_idx_end+1] * volume_lat, axis=1) / np.sum(volume_lat, axis=1), depth,  label='N² (1-100)', color='blue')
plt.plot(np.nansum(N2_lat_500[:, lat_idx:lat_idx_end+1] * volume_lat, axis=1) / np.sum(volume_lat, axis=1), depth, label='N² (500-600)', color='red')
#plt.plot(np.mean(dN2, axis=1), depth, label='ΔN²', color='black')
plt.xlabel('Depth [m]', fontsize=12)    
plt.ylim(2000, 0)

plt.figure(figsize=(8, 6))
plt.plot(N2_som_1_WGKP, depth,  label='N² (1-100)', color='blue')
plt.plot(N2_som_2_WGKP, depth, label='N² (500-600)', color='red')
#plt.plot(np.mean(dN2, axis=1), depth, label='ΔN²', color='black')
plt.xlabel('Depth [m]', fontsize=12) 
plt.title('Over convection region (35$^\circ$W - 80$^\circ$W, 80 - 50$^\circ$S)', fontsize=14)   
plt.ylim(2000, 0)

# %%

plt.figure(figsize=(8, 6))
plt.contourf(lat, depth, N2_lat_500 - N2_lat_1, levels= np.linspace(-2e-6, 2e-6, 21), extend = 'both', cmap = 'RdBu_r')
cbar	= colorbar(CS2, ticks = np.arange(-2e-6, 2e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
plt.title('Buoyancy Frequency difference (500-600) vs (1-100)', fontsize=14)
plt.xlim(-80,-30)
plt.ylim(500, 0)
plt.contour(lat, depth, dens_1, levels = [1027.1], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [1027.1], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [1027.5], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [1027.5], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [1027.7], colors = 'black')
plt.contour(lat, depth, dens_500, levels = [1027.7], linestyles = '--', colors='black')
plt.ylabel('Depth [m]', fontsize=12)
plt.xlabel('Latitude [$^\circ$S]', fontsize=12)
plt.vlines(x=-80, ymin=2000, ymax=0, color='grey')
plt.vlines(x=-50, ymin=2000, ymax=0, color='grey')

# %%
