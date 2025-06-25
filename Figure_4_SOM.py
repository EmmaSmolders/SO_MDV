#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 12:24:48 2025

@author: 6008399

Figure 4 SOM paper

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

#%% Functions

def compute_N2(dens_before, dens_after, depth):
    # Centered difference
    density_before = dens_before
    drho_dz_before = np.zeros_like(density_before)
    drho_dz_before[1:-1] = (density_before[2:] - density_before[:-2]) / (depth[2:] - depth[:-2])

    density_after = dens_after
    drho_dz_after = np.zeros_like(density_after)
    drho_dz_after[1:-1] = (density_after[2:] - density_after[:-2]) / (depth[2:] - depth[:-2])

    # Forward/backward difference at boundaries
    drho_dz_before[0] = (density_before[1] - density_before[0]) / (depth[1] - depth[0])
    drho_dz_before[-1] = (density_before[-1] - density_before[-2]) / (depth[-1] - depth[-2])

    drho_dz_after[0] = (density_after[1] - density_after[0]) / (depth[1] - depth[0])
    drho_dz_after[-1] = (density_after[-1] - density_after[-2]) / (depth[-1] - depth[-2])

    # Compute N² = −(g / ρ) (dρ/dz)
    g = 9.81  # m/s^2
    N2_before = - (g / density_before) * drho_dz_before
    N2_after = - (g / density_after) * drho_dz_after

    return N2_before, N2_after

#%% Mixed layer depth WGKP

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_year_1-600_area_averaged_WGKP.nc','r')

time_temp_wgkp = fh.variables['time'][:] #Time SOM cycle 1
depth_temp_wgkp= fh.variables['depth'][:] 
temp_wgkp = fh.variables['TEMP'][:]

fh.close()

fh = netcdf.Dataset(directory + 'PD_year_2-600_area_averaged_WGKP.nc','r')

time_PD_wgkp = fh.variables['time'][:] #Time SOM cycle 1
depth_PD_wgkp= fh.variables['depth'][:] 
PD_wgkp = fh.variables['PD'][:]

fh.close()

#%% WGKP densities

fh = netcdf.Dataset(directory + 'PD_year_63-114_area_averaged_WGKP.nc','r')

depth       = fh.variables['depth'][:]  #depth
time        = fh.variables['time'][:]
dens_1_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'PD_year_324-378_area_averaged_WGKP.nc','r')

depth       = fh.variables['depth'][:]  #depth
time        = fh.variables['time'][:]
dens_2_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'PD_year_500-600_area_averaged_WGKP.nc','r')

depth       = fh.variables['depth'][:]  #depth
time        = fh.variables['time'][:]
dens_3_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()

#%% Compute N2

N2_before_wgkp_1_2, N2_after_wgkp_1_2 = compute_N2(np.mean(dens_1_wgkp, axis=0), np.mean(dens_2_wgkp, axis=0), depth)
N2_before_wgkp_1_3, N2_after_wgkp_1_3 = compute_N2(np.mean(dens_1_wgkp, axis=0), np.mean(dens_3_wgkp, axis=0), depth)

#%%

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))

CS = ax1.contourf(time_temp_wgkp, depth_temp_wgkp, (temp_wgkp - np.mean(temp_wgkp[0:50], axis=0)).transpose(), levels=np.linspace(-0.4, .4, 21), extend='both', cmap = 'RdBu_r')
cbar	= colorbar(CS)
cbar.set_label('Temperature anomaly [$^\circ$C]', fontsize = 11)
ax1.plot(time_mld, MLD_max, color='black', label='maximum MLD')
ax1.set_title('a) MLD maximum and temperature anomaly', fontsize=14)
ax1.set_ylim(5000, 1)
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Time [model years]', fontsize=12)
ax1.legend(loc=3, fontsize=11)

ax2.plot(-N2_after_wgkp_1_2 - -N2_before_wgkp_1_2, depth, color='blue', label='SOM2 vs SOM1')
ax2.plot(-N2_after_wgkp_1_3 - -N2_before_wgkp_1_3, depth, color='red', label='SOM3 vs SOM1')
ax2.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax2.set_title('b) Area-averaged N$^2$ difference', fontsize=14)
ax2.set_xlabel('N$^2$ difference [s$^{-1}$]', fontsize=12)
#ax2.set_ylabel('Depth [m]')
ax2.legend(loc=4, fontsize=11)
ax2.set_ylim(5000, 1)
ax2.set_xlim(-1.e-6, 4e-6)

fig.savefig(directory_figures +'MLD_WGKP_max_1-600_N2_before_after.pdf')

#%%

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))

CS = ax1.contourf(time_PD_wgkp, depth_PD_wgkp, (PD_wgkp - np.mean(PD_wgkp[50:150], axis=0)).transpose(), levels=np.linspace(-0.06, .06, 21), extend='both', cmap = 'RdBu_r')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-0.06, -0.04, -0.02,  0,  0.02,  0.04, 0.06])  # Set specific ticks
cbar.set_label('Potential density anomaly [kg/m$^3$]', fontsize = 11)
ax1.plot(time_mld, MLD_max, color='black', label='maximum MLD')
ax1.set_title('a) MLD maximum and potential density anomaly', fontsize=14)
ax1.set_ylim(5000, 1)
ax1.set_xlim(35,600)
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Time [model years]', fontsize=12)
ax1.legend(loc=3, fontsize=11)

ax2.plot(-N2_after_wgkp_1_2 - -N2_before_wgkp_1_2, depth, color='blue', label='SOM2 vs SOM1')
ax2.plot(-N2_after_wgkp_1_3 - -N2_before_wgkp_1_3, depth, color='red', label='SOM3 vs SOM1')
ax2.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
ax2.set_title('b) Area-averaged N$^2$ difference', fontsize=14)
ax2.set_xlabel('N$^2$ difference [s$^{-1}$]', fontsize=12)
#ax2.set_ylabel('Depth [m]')
ax2.legend(loc=4, fontsize=11)
ax2.set_ylim(5000, 1)
ax2.set_xlim(-1.e-6, 4e-6)

fig.savefig(directory_figures +'PD_MLD_WGKP_max_1-600_N2_before_after.pdf')

#%% Model year 28 to 35 contains in-situ density, insead of potential density. Discard this.

plt.figure()
plt.contourf(time_PD_wgkp, depth_PD_wgkp, PD_wgkp.transpose())
plt.colorbar()


