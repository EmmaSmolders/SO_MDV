#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 12 15:40:51 2025

@author: 6008399

APE for different sectors of the Southern Ocean

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

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-600_window_5_volume_integrated_SO30_TEST_rhoref_0-50_TEST_goedegsw.nc','r')

time_1_600 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_SO30 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_350-600_window_5_volume_integrated_SO30_TEST_rhoref_0-50_TEST_goedegsw.nc','r')

time_350_600 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_SO30_350_600 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_63-114_window_5_volume_integrated_SO30_TEST_rhoref_63-114_TEST_goedegsw.nc','r')

time_som_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_SO30_som_1 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_410-480_window_5_volume_integrated_SO30_TEST_rhoref_410-480_TEST_goedegsw.nc','r')

time_som_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_SO30_som_2 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_500-600_window_5_volume_integrated_SO30_TEST_rhoref_500-600_TEST_goedegsw.nc','r')

time_som_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_SO30_som_3 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-200_window_5_volume_integrated_Pacific_rhoref_0-50.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Pacific_1 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_190-400_window_5_volume_integrated_Pacific_rhoref_0-50.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Pacific_2 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_390-600_window_5_volume_integrated_Pacific_rhoref_0-50.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Pacific_3 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-200_window_5_volume_integrated_entire_Pacific_rhoref_0-50.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_entire_Pacific_1 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_190-400_window_5_volume_integrated_entire_Pacific_rhoref_0-50.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_entire_Pacific_2 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_390-600_window_5_volume_integrated_entire_Pacific_rhoref_0-50.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_entire_Pacific_3 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-200_window_5_volume_integrated_WGKP_rhoref_0-50.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_wgkp_1 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_190-400_window_5_volume_integrated_WGKP_rhoref_0-50.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_wgkp_2 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_390-600_window_5_volume_integrated_WGKP_rhoref_0-50.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_wgkp_3 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-200_window_5_volume_integrated_Atlantic_rhoref_0-50.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Atlantic_1 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_190-400_window_5_volume_integrated_Atlantic_rhoref_0-50.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Atlantic_2 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_390-600_window_5_volume_integrated_Atlantic_rhoref_0-50.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Atlantic_3 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-200_window_5_volume_integrated_Indian_rhoref_0-50.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Indian_1 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_190-400_window_5_volume_integrated_Indian_rhoref_0-50.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Indian_2 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_390-600_window_5_volume_integrated_Indian_rhoref_0-50.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_Indian_3 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data

#%%

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-600_window_5_volume_integrated_AU_rhoref_0-50.nc','r')

time_all = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_AU = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-600_window_5_volume_integrated_NZ_rhoref_0-50.nc','r')

time_all = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_NZ = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/APE_year_1-600_window_5_volume_integrated_PA_rhoref_0-50.nc','r')

time_all = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_PA = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%% Total energetics

fig, (ax1, ax3, ax4, ax5) = plt.subplots(4, 1, figsize=(8, 10))

ax1.plot(time_som_1, APE_SO30_som_1/1e18, label='$\\rho_{ref}$ SOM', color='cyan')
ax1.plot(time_som_2, APE_SO30_som_2/1e18, color='cyan')
ax1.plot(time_som_3, APE_SO30_som_3/1e18, color='cyan')
ax1.plot(time_1_600, APE_SO30/1e18, label='$\\rho_{ref}$ 0-50', color='royalblue')
ax1.plot(time_350_600, APE_SO30_350_600/1e18, color='royalblue')
ax1.set_title('a) SO30')
ax1.set_ylabel('P [EJ]')
ax1.tick_params(labelcolor='black')
ax1.grid(True)
ax1.set_xlim(0,600)
ax1.legend()

#ax2.plot(time_1, APE_wgkp_1/1e18, color='red', label='WGKP')
#ax2.plot(time_2, APE_wgkp_2/1e18, color='red')
#ax2.plot(time_3, APE_wgkp_3/1e18, color='red')
#ax2.set_title('b) WGKP (35$^\circ$W - 80$^\circ$E, 80 - 50$^\circ$S)')#' and depth averaged SOM (0-1000m)')
#ax2.set_ylabel('P [EJ]', color='black')
#ax2.tick_params(axis='y', labelcolor='black')
#ax2.grid(True)

ax3.plot(time_1, APE_Atlantic_1/1e18, color='black', label='Atlantic')
ax3.plot(time_2, APE_Atlantic_2/1e18, color='black')
ax3.plot(time_3, APE_Atlantic_3/1e18, color='black')
ax3.set_xlim(0,600)
ax3.set_title('c) Atlantic sector (65$^\circ$W - 25$^\circ$E)')#' and depth averaged SOM (0-1000m)')
ax3.set_ylabel('P [EJ]', color='black')
ax3.tick_params(axis='y', labelcolor='black')
ax3.grid(True)

#ax3.plot(time_1, APE_Pacific_1/1e18, color='green', label='Pacific')
#ax3.plot(time_2, APE_Pacific_2/1e18, color='green')
#ax3.plot(time_3, APE_Pacific_3/1e18, color='green')

ax4.plot(time_1, APE_entire_Pacific_1/1e18, color='green', label='entire Pacific')
ax4.plot(time_2, APE_entire_Pacific_2/1e18, color='green')
ax4.plot(time_3, APE_entire_Pacific_3/1e18, color='green')
ax4.set_title('d) Pacific sector (150$^\circ$E - 65$^\circ$W)')
ax4.set_ylabel('P [EJ]', color='black')
ax4.set_xlim(0,600)
ax4.grid(True)
#ax4.legend()

ax5.plot(time_1, APE_Indian_1/1e18, color='magenta', label='Indian')
ax5.plot(time_2, APE_Indian_2/1e18, color='magenta')
ax5.plot(time_3, APE_Indian_3/1e18, color='magenta')
ax5.set_title('e) Indian sector (25$^\circ$E - 150$^\circ$E)')
ax5.set_ylabel('P [EJ]', color='black')
ax5.set_xlim(0,600)
ax5.grid(True)
#ax4.legend()


plt.suptitle('Energetics SO30 and its basins (90 - 30$^\circ$S)', fontsize=16)
plt.tight_layout()
plt.savefig(directory_figures +'APE_SO30_basins.pdf')
plt.show()

#%%

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 10))

ax1.plot(time_all, APE_NZ/1e18, color='royalblue')
ax1.set_title('a) NZ convective region')
ax1.set_ylabel('P [EJ]')
ax1.tick_params(labelcolor='black')
ax1.grid(True)
ax1.set_xlim(0,600)
#ax1.legend()

ax2.plot(time_1, APE_wgkp_1/1e18, color='red', label='WGKP')
ax2.plot(time_2, APE_wgkp_2/1e18, color='red')
ax2.plot(time_3, APE_wgkp_3/1e18, color='red')
ax2.set_title('b) WGKP convective region')#' and depth averaged SOM (0-1000m)')
ax2.set_ylabel('P [EJ]', color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)

ax3.plot(time_all, APE_AU/1e18, color='black', label='Atlantic')
ax3.set_xlim(0,600)
ax3.set_title('c) AU convective region')#' and depth averaged SOM (0-1000m)')
ax3.set_ylabel('P [EJ]', color='black')
ax3.tick_params(axis='y', labelcolor='black')
ax3.grid(True)

#ax3.plot(time_1, APE_Pacific_1/1e18, color='green', label='Pacific')
#ax3.plot(time_2, APE_Pacific_2/1e18, color='green')
#ax3.plot(time_3, APE_Pacific_3/1e18, color='green')

ax4.plot(time_all, APE_PA/1e18, color='green', label='entire Pacific')
ax4.set_title('d) PA convective region')
ax4.set_ylabel('P [EJ]', color='black')
ax4.set_xlim(0,600)
ax4.grid(True)
#ax4.legend()

plt.suptitle('APE of convective regions', fontsize=16)
plt.tight_layout()
plt.savefig(directory_figures +'APE_convective regions.pdf')
plt.show()
