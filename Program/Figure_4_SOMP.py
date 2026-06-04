#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 27 22:12:41 2025

@author: 6008399

SOMP Figure 4

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
directory = r'/Users/6008399/Documents/PhD/HR_POP/netcdf/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%% Functions

def compute_N2_from_profile(density, depth, rho0=1027, g=9.81):
    """
    Compute vertical density gradient (drho/dz) and buoyancy frequency N²
    from a 1D density profile and depth array.

    Compute N² exactly following the user's original code,
    assuming depth increases downward.
    """

    nz = len(depth)
    n0 = np.zeros(nz)

    for i in range(nz):
        if i == 0:
            # surface: (rho0 - rho1)/(z1 - z0)
            n0[i] = (density[i] - density[i+1]) / (depth[i+1] - depth[i])

        elif i == nz - 1:
            # bottom: (rho_{n-2} - rho_{n-1})/(z_{n-1}-z_{n-2})
            n0[i] = (density[i-1] - density[i]) / (depth[i] - depth[i-1])

        else:
            # centered: (rho[i-1] - rho[i+1]) / (z[i+1] - z[i-1])
            n0[i] = (density[i-1] - density[i+1]) / (depth[i+1] - depth[i-1])

    # Your formula: N² = -(g / rho0) * n0
    N2 = - (g / rho0) * n0

    return n0, N2

#%% Mixed layer depth NZ

fh = netcdf.Dataset(directory + 'Ocean/MLD_year_1-600_NZ.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_NZ = fh.variables['MLD_max'][:] 
MLD = fh.variables['MLD'][:]

fh.close()

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_wgkp= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'PD_year_2-600_area_averaged_WGKP.nc','r')

time_PD_wgkp = fh.variables['time'][:] #Time SOM cycle 1
depth_PD_wgkp= fh.variables['depth'][:] 
PD_WGKP = fh.variables['PD'][:]

fh.close()


fh = netcdf.Dataset(directory + 'PD_year_2-600_area_averaged_AU.nc','r')

depth       = fh.variables['depth'][:]  #depth
time_PD_AU = fh.variables['time'][:] #Time SOM cycle 1
depth_PD_AU= fh.variables['depth'][:] 
PD_AU = fh.variables['PD'][:]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/PD_year_2-600_area_averaged_NZ_70_60S_150_190E.nc','r')

depth       = fh.variables['depth'][:]  #depth
time_PD_NZ = fh.variables['time'][:] #Time SOM cycle 1
depth_PD_NZ= fh.variables['depth'][:] 
PD_NZ = fh.variables['PD'][:]

fh.close()

#%% NZ densities

time_1              = time_PD_NZ[62:113]
dens_1_NZ         = PD_NZ[62:113,:]     #potential density

time_2              = time_PD_NZ[410:480]
dens_2_NZ         = PD_NZ[410:480,:]     #potential density

time_3              = time_PD_NZ[499::]
dens_3_NZ         = PD_NZ[499::,:]     #potential density

dens_1_WGKP         = PD_WGKP[62:113,:]     #potential density
dens_2_WGKP         = PD_WGKP[410:480,:]     #potential density
dens_3_WGKP         = PD_WGKP[499::,:]     #potential density

plt.figure()
plt.plot(np.mean(dens_3_NZ, axis=0), depth)

#%%Test
n0 = np.zeros(len(depth))
for depth_i in range(len(depth)):
	if depth_i == 0:
		n0[depth_i] = (np.mean(dens_1_NZ[:,depth_i], axis=0) - np.mean(dens_1_NZ[:,depth_i+1], axis=0))/(depth[1] - depth[0])
	elif depth_i == len(depth) - 1:
		n0[depth_i] = (np.mean(dens_1_NZ[:,depth_i - 1], axis=0) - np.mean(dens_1_NZ[:,depth_i], axis=0))/(depth[-1] - depth[-2])
	else:
		n0[depth_i] = (np.mean(dens_1_NZ[:,depth_i-1], axis=0) - np.mean(dens_1_NZ[:,depth_i+1], axis=0))/(np.abs(depth[depth_i-1] - depth[depth_i+1]))

plt.figure()
plt.plot(n0, depth)

rho0 = 1027
g = 9.81

N2 = - g/rho0 * n0

#Something weird happens in the lower layers (PD profiles are also reaally weird here. Upper 3000m is correct though)
plt.figure()
plt.plot(N2, depth)

n0_som_1_NZ, N2_som_1_NZ = compute_N2_from_profile(np.mean(dens_1_NZ, axis=0), depth)
n0_som_2_NZ, N2_som_2_NZ = compute_N2_from_profile(np.mean(dens_2_NZ, axis=0), depth)
n0_som_3_NZ, N2_som_3_NZ = compute_N2_from_profile(np.mean(dens_3_NZ, axis=0), depth)

n0_som_1_WGKP, N2_som_1_WGKP = compute_N2_from_profile(np.mean(dens_1_WGKP, axis=0), depth)
n0_som_2_WGKP, N2_som_2_WGKP = compute_N2_from_profile(np.mean(dens_2_WGKP, axis=0), depth)
n0_som_3_WGKP, N2_som_3_WGKP = compute_N2_from_profile(np.mean(dens_3_WGKP, axis=0), depth)

plt.figure()
plt.plot(N2_som_1_NZ, depth)

plt.figure()
plt.plot(N2_som_1_WGKP, depth)

#%%

def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

window = 5

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS = ax1.contourf(time_PD_wgkp, depth_PD_wgkp, (PD_WGKP - np.mean(PD_WGKP[50:150], axis=0)).transpose(), levels=np.linspace(-0.06, .06, 21), extend='both', cmap = 'RdBu_r')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-0.06, -0.04, -0.02,  0,  0.02,  0.04, 0.06])  # Set specific ticks
cbar.set_label('Potential density anomaly [kg/m$^3$]', fontsize = 11)
ax1.plot(time_mld[window//2:-window//2+1], Moving_average(MLD_max_wgkp, window), color='black', label='maximum MLD')
ax1.set_title('a) MLD maximum and PD anomaly (WGKP)', fontsize=14)
#ax1.plot(time_mld, np.mean(MLD, axis=(1,2)))
ax1.set_ylim(2000, 1)
ax1.set_xlim(35,600)
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Time [model years]', fontsize=12)
ax1.legend(loc=4, fontsize=11)

ax2.plot(N2_som_2_WGKP - N2_som_1_WGKP, depth, color='blue', label='SOM2 vs SOM1')
ax2.plot(N2_som_3_WGKP - N2_som_1_WGKP, depth, color='red', label='SOM3 vs SOM1')
ax2.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax2.set_title('b) Area-averaged N$^2$ difference (WGKP)', fontsize=14)
ax2.set_xlabel('N$^2$ difference [s$^{-1}$]', fontsize=12)
#ax2.set_ylabel('Depth [m]')
ax2.legend(loc=4, fontsize=11)
ax2.set_ylim(5000, 1)
ax2.set_ylim(2000, 1)
ax2.set_xlim(-4.e-6, 10e-6)

CS = ax3.contourf(time_PD_NZ, depth_PD_NZ, (PD_NZ - np.mean(PD_NZ[50:150], axis=0)).transpose(), levels=np.linspace(-0.06, .06, 21), extend='both', cmap = 'RdBu_r')
cbar = fig.colorbar(CS, ax=ax3, ticks=[-0.06, -0.04, -0.02,  0,  0.02,  0.04, 0.06])  # Set specific ticks
cbar.set_label('Potential density anomaly [kg/m$^3$]', fontsize = 11)
ax3.plot(time_mld[window//2:-window//2+1], Moving_average(MLD_max_NZ, window), color='black', label='maximum MLD')
ax3.set_title('c) MLD maximum and PD anomaly (Pacific)', fontsize=14)
#ax1.plot(time_mld, np.mean(MLD, axis=(1,2)))
ax3.set_ylim(2000, 1)
ax3.set_xlim(35,600)
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Time [model years]', fontsize=12)
ax3.legend(loc=4, fontsize=11)

ax4.plot(N2_som_2_NZ - N2_som_1_NZ, depth, color='blue', label='SOM2 vs SOM1')
ax4.plot(N2_som_3_NZ - N2_som_1_NZ, depth, color='red', label='SOM3 vs SOM1')
ax4.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax4.set_title('d) Area-averaged N$^2$ difference (Pacific)', fontsize=14)
ax4.set_xlabel('N$^2$ difference [s$^{-1}$]', fontsize=12)
#ax2.set_ylabel('Depth [m]')
ax4.legend(loc=4, fontsize=11)
ax4.set_ylim(5000, 1)
ax4.set_ylim(2000, 1)
ax4.set_xlim(-4.e-6, 10e-6)

#plt.suptitle('Convection region Pacific', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Figure_6_SOM_OS.pdf')

#%% Model year 28 to 35 contains in-situ density, insead of potential density. Discard this.

plt.figure()
plt.contourf(time_PD_NZ, depth_PD_NZ, PD_NZ.transpose())
plt.colorbar()

#%%

plt.figure()
plt.plot(Moving_average(MLD_max_NZ, window))
plt.plot(Moving_average(np.mean(MLD, axis=(1,2)), window))


