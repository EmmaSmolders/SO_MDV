#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 15:48:05 2025

@author: 6008399

Figure 1 SOM paper

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

#%% Read in data 

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh1 = netcdf.Dataset(directory + 'AMOC_transport_depth_0-1000m.nc','r')

amoc_strength   = fh1.variables['Transport'][:]     #Volume transport [Sv]

fh1.close()

fh2 = netcdf.Dataset(directory + 'DSL_Hovmoller.nc','r')

lat = fh2.variables['lat'][:] #Latitude [deg N]
DSL = fh2.variables['DSL'][:] #Dynamic sea level [cm]

fh2.close()

fh3 = netcdf.Dataset(directory + 'TEMP_SALT_Hovmoller_depth_300-700m.nc','r')

temp = fh3.variables['TEMP'][:] #Temperature [degC]
salt = fh3.variables['SALT'][:] #Salinity [g/kg]

fh3.close()

fh = netcdf.Dataset(directory + 'Drake_Passage_transport.nc','r')

time                = fh.variables['time'][:]         #time [model years]
acc_transport           = fh.variables['Transport'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_depth_0-1000m.nc','r')

time                = fh.variables['time'][:]         #time [model years]
som_star           = fh.variables['SOM'][:]    #Transport [Sv]

fh.close()

#%% Plot results

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8))

ax1.plot(time*0.0003, amoc_strength, label='AMOC Strength', color='blue')
ax1.set_title('a) AMOC strength at 26$^\circ$N')
#ax1.set_xlabel('Time [model years]')
ax1.set_ylabel('Volume transport [Sv]')
ax1.tick_params(labelcolor='black')
ax1.grid(True)
ax1.set_xlim(0,0.18)

ax2.plot(time*0.0003, som, label='SOM', color='black')
ax2.set_title('b) Southern Ocean Mode (SOM)')#' and depth averaged SOM (0-1000m)')
#ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('SOM [$^\circ$C]', color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(time*0.0003, som_star, label='SOM$^*$', color='red')
ax2_twin.set_ylabel('SOM$^*$[$^\circ$C]', color='red')
ax2_twin.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(0,0.18)

ax3.plot(time*0.0003, acc_transport, color='green')
ax3.set_title('c) Drake Passage transport')
ax3.set_xlabel('Freshwater flux forcing F$_H$ [Sv]')
ax3.set_ylabel('Volume transport [Sv]')
ax3.set_ylim(94,125)
ax3.set_xlim(0,0.18)
ax3.grid(True)

ax4 = fig.add_axes([0.7, 0.205, 0.302, 0.150], projection=ccrs.SouthPolarStereo())
ax4.set_extent([-50, 60, -90, -33], crs=ccrs.PlateCarree())

ax4.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
ax4.add_feature(cfeature.OCEAN, zorder=0, edgecolor='black')

gl = ax4.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = False

#Define region for som regions
lon_min, lon_max = -50, 0
lat_min, lat_max = -50, -35

lon_min_star, lon_max_star = -30, 20
lat_min_star, lat_max_star = -55, -40

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

ax4.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_1, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_2, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_star	= np.arange(lon_min_star, lon_max_star + 0.1, 0.1)
y_1_star	= np.zeros(len(x_1)) + lat_min_star
y_2_star	= np.zeros(len(x_1)) + lat_max_star

y_3_star = np.arange(lat_min_star, lat_max_star + 0.1, 0.1)
x_2_star = np.zeros(len(y_3)) + lon_min_star
x_3_star = np.zeros(len(y_3)) + lon_max_star

ax4.plot(x_1_star, y_1_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_1_star, y_2_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_2_star, y_3_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_3_star, y_3_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_acc	= np.arange(acc_lon_min, acc_lon_max + 0.1, 0.1)
y_1_acc	= np.zeros(len(x_1_acc)) + acc_lat_min
y_2_acc	= np.zeros(len(x_1_acc)) + acc_lat_max

y_3_acc = np.arange(acc_lat_min, acc_lat_max + 0.1, 0.1)
x_2_acc = np.zeros(len(y_3_acc)) + -66
x_3_acc = np.zeros(len(y_3_acc)) + acc_lon_max

ax4.plot(x_1_acc, y_1_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_1_acc, y_2_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_2_acc, y_3_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp

ax4.plot(x_1_wgkp, y_1_wgkp, '-', color='blue', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_1_wgkp, y_2_wgkp, '-', color='blue', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_2_wgkp, y_3_wgkp, '-', color='blue', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_3_wgkp, y_3_wgkp, '-', color='blue', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

plt.tight_layout()
plt.savefig(directory_figures +'SOM_AMOC_transport_HR_pop.pdf')
plt.show()
