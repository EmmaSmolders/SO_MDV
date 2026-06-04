#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 12:15:54 2025

@author: 6008399

Figure S1 of SOM paper

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
directory = r'/Users/6008399/Documents/PhD/HR_POP/Zenodo/Data_Final/'
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

fh3 = netcdf.Dataset(directory + 'TEMP_SALT_Hovmoller_depth_300-700m_Pacific_170_290_E.nc','r')

temp_somp = fh3.variables['TEMP'][:] #Temperature [degC]
salt_somp = fh3.variables['SALT'][:] #Salinity [g/kg]

fh3.close()

fh = netcdf.Dataset(directory + 'Drake_Passage_transport.nc','r')

time                = fh.variables['time'][:]         #time [model years]
acc_transport           = fh.variables['Transport'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_depth_0-1000m.nc','r')

time                = fh.variables['time'][:]         #time [model years]
som_star           = fh.variables['SOM'][:]    #Transport [Sv]

fh.close()

#%%

DSL_anomaly     = np.zeros((len(time), len(lat)))
temp_anomaly    = np.zeros((len(time), len(lat)))
salt_anomaly    = np.zeros((len(time), len(lat)))
temp_anomaly_somp    = np.zeros((len(time), len(lat)))
salt_anomaly_somp    = np.zeros((len(time), len(lat)))

for lat_i in range(len(lat)):    
    DSL_anomaly[:, lat_i] = DSL[:, lat_i] - np.mean(DSL[0:100,lat_i])
    temp_anomaly[:, lat_i] = temp[:, lat_i] - np.mean(temp[0:100,lat_i])
    salt_anomaly[:, lat_i] = salt[:, lat_i] - np.mean(salt[0:100,lat_i])
    temp_anomaly_somp[:, lat_i] = temp_somp[:, lat_i] - np.mean(temp_somp[0:100,lat_i])
    salt_anomaly_somp[:, lat_i] = salt_somp[:, lat_i] - np.mean(salt_somp[0:100,lat_i])

plt.figure()
plt.contourf(time, lat, DSL_anomaly.transpose(), cmap='seismic', levels=np.linspace(-21.1, 21.1,21), extend='both')
plt.ylabel('Latitude [$^\circ$N]')
plt.colorbar(label='Sea surface height anomaly [cm]')
plt.title('SSH anomalies, HR-POP')
plt.xlabel('Time [model years]')
plt.tight_layout()
plt.savefig(directory_figures +'SSH_anomalies_HR_pop.pdf')
plt.show()

#%%

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(13, 5))

CS = ax1.contourf(time, lat, temp_anomaly.transpose(), cmap='seismic', levels=np.linspace(-5., 5., 21), extend='both')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])  # Set specific ticks
cbar.set_label(label = 'Temperature anomaly [$^\circ$C]', fontsize=14)
ax1.set_ylabel('Latitude [$^\circ$N]', fontsize=16)
ax1.set_xlabel('Time [model years]', fontsize=16)
ax1.set_title('a) Temperature anomalies', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
cbar.ax.tick_params(labelsize=13)

CS1 = ax2.contourf(time, lat, salt_anomaly.transpose(), cmap='BrBG_r', levels=np.linspace(-0.5, 0.5, 21), extend='both')
cbar1 = fig.colorbar(CS1, ax=ax2, ticks=[-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5])  # Set specific ticks
cbar1.ax.tick_params(labelsize=13)
cbar1.set_label(label = 'Salinity anomaly [g/kg]', fontsize=14)
#ax2.set_ylabel('Latitude [$^\circ$N]')
ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_title('b) Salinity anomalies', fontsize=17)
ax2.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig(directory_figures +'SALT_TEMP_SOM_anomalies_HR_pop.pdf')
plt.show()

#%%

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(13, 5))

CS = ax1.contourf(time[0:200], lat, temp_anomaly[0:200].transpose(), cmap='seismic', levels=np.linspace(-2., 2., 21), extend='both')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2])  # Set specific ticks
cbar.set_label(label = 'Temperature anomaly [$^\circ$C]', fontsize=14)
ax1.set_ylabel('Latitude [$^\circ$N]', fontsize=16)
ax1.set_xlabel('Time [model years]', fontsize=16)
ax1.set_title('a) Temperature anomalies', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
cbar.ax.tick_params(labelsize=13)

CS1 = ax2.contourf(time[0:200], lat, salt_anomaly[0:200].transpose(), cmap='BrBG_r', levels=np.linspace(-0.1, 0.1, 21), extend='both')
cbar1 = fig.colorbar(CS1, ax=ax2, ticks=[-0.1, 0, 0.1])  # Set specific ticks
cbar1.ax.tick_params(labelsize=13)
cbar1.set_label(label = 'Salinity anomaly [g/kg]', fontsize=14)
#ax2.set_ylabel('Latitude [$^\circ$N]')
ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_title('b) Salinity anomalies', fontsize=17)
ax2.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig(directory_figures +'SALT_TEMP_SOM_anomalies_HR_pop_200years.pdf')
plt.show()
# %%

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(13, 5))

CS = ax1.contourf(time, lat, temp_anomaly_somp.transpose(), cmap='seismic', levels=np.linspace(-1., 1., 21), extend='both')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1])  # Set specific ticks
cbar.set_label(label = 'Temperature anomaly [$^\circ$C]', fontsize=14)
ax1.set_ylabel('Latitude [$^\circ$N]', fontsize=16)
ax1.set_xlabel('Time [model years]', fontsize=16)
ax1.set_title('a) Temperature anomalies', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
cbar.ax.tick_params(labelsize=13)

CS1 = ax2.contourf(time, lat, salt_anomaly_somp.transpose(), cmap='BrBG_r', levels=np.linspace(-0.1, 0.1, 21), extend='both')
cbar1 = fig.colorbar(CS1, ax=ax2, ticks=[-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1])  # Set specific ticks
cbar1.ax.tick_params(labelsize=13)
cbar1.set_label(label = 'Salinity anomaly [g/kg]', fontsize=14)
#ax2.set_ylabel('Latitude [$^\circ$N]')
ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_title('b) Salinity anomalies', fontsize=17)
ax2.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig(directory_figures +'SALT_TEMP_SOMP_anomalies_HR_pop.pdf')
plt.show()
# %% Longitudinal hovmoller diagrams 

fh3 = netcdf.Dataset(directory + 'TEMP_SALT_Hovmoller_depth_0-5000m_lat_-58--40N.nc','r')

temp_lon_58_40S = fh3.variables['TEMP'][:] #Temperature [degC]
salt_lon_58_40S = fh3.variables['SALT'][:] #Salinity [g/kg]
lon = fh3.variables['lon'][:] #Longitude [deg E] 
time = fh3.variables['time'][:] #Time [model years]   

fh3.close()

fh3 = netcdf.Dataset(directory + 'TEMP_SALT_Hovmoller_depth_0-5000m_lat_-70--58N.nc','r')

temp_lon_70_58S = fh3.variables['TEMP'][:] #Temperature [degC]
salt_lon_70_58S = fh3.variables['SALT'][:] #Salinity [g/kg]

fh3.close()

#%%

plt.figure()
plt.contourf(temp_lon_58_40S.transpose())
plt.ylabel('Longitude [$^\circ$E]')
plt.colorbar(label='Temperature [degC]')
plt.title('Temperature, 58-40S, HR-POP')
plt.xlabel('Time [model years]')
plt.tight_layout()
#plt.savefig(directory_figures +'TEMP_lon_58_40S_HR_pop.pdf')
plt.show()

plt.figure()
plt.contourf(time, lon, temp_lon_70_58S.transpose())
plt.ylabel('Longitude [$^\circ$E]')
plt.colorbar(label='Temperature [degC]')
plt.title('Temperature, 70-58S, HR-POP')
plt.xlabel('Time [model years]')
plt.tight_layout()
#plt.savefig(directory_figures +'TEMP_lon_58_40S_HR_pop.pdf')
plt.show()

#%%

temp_lon_58_40S_anomaly    = np.zeros((len(time), len(lon)))
salt_lon_58_40S_anomaly    = np.zeros((len(time), len(lon)))
temp_lon_70_58S_anomaly    = np.zeros((len(time), len(lon)))
salt_lon_70_58S_anomaly    = np.zeros((len(time), len(lon)))

for lon_i in range(len(lon)):    
    temp_lon_58_40S_anomaly[:, lon_i] = temp_lon_58_40S[:, lon_i] - np.mean(temp_lon_58_40S[0:100,lon_i])
    salt_lon_58_40S_anomaly[:, lon_i] = salt_lon_58_40S[:, lon_i] - np.mean(salt_lon_58_40S[0:100,lon_i])
    temp_lon_70_58S_anomaly[:, lon_i] = temp_lon_70_58S[:, lon_i] - np.mean(temp_lon_70_58S[0:100,lon_i])
    salt_lon_70_58S_anomaly[:, lon_i] = salt_lon_70_58S[:, lon_i] - np.mean(salt_lon_70_58S[0:100,lon_i])

# %%

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(13, 5))

CS = ax1.contourf(temp_lon_58_40S_anomaly, cmap='seismic', levels=np.linspace(-0.05, 0.05, 21), extend='both')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1])  # Set specific ticks
cbar.set_label(label = 'Temperature anomaly [$^\circ$C]', fontsize=14)
ax1.set_ylabel('Latitude [$^\circ$N]', fontsize=16)
ax1.set_xlabel('Time [model years]', fontsize=16)
ax1.set_title('a) Temperature anomalies', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
cbar.ax.tick_params(labelsize=13)

CS1 = ax2.contourf(salt_lon_58_40S_anomaly, cmap='BrBG_r', levels=np.linspace(-0.001, 0.001, 21), extend='both')
cbar1 = fig.colorbar(CS1, ax=ax2, ticks=[-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1])  # Set specific ticks
cbar1.ax.tick_params(labelsize=13)
cbar1.set_label(label = 'Salinity anomaly [g/kg]', fontsize=14)
#ax2.set_ylabel('Latitude [$^\circ$N]')
ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_title('b) Salinity anomalies', fontsize=17)
ax2.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig(directory_figures +'SALT_TEMP_anomalies_lon_58_40S_HR_pop.pdf')
plt.show()
# %%
fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(13, 5))

CS = ax1.contourf(temp_lon_70_58S_anomaly, cmap='seismic', levels=np.linspace(-0.05, 0.05, 21), extend='both')
cbar = fig.colorbar(CS, ax=ax1, ticks=[-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1])  # Set specific ticks
cbar.set_label(label = 'Temperature anomaly [$^\circ$C]', fontsize=14)
ax1.set_xlabel('Longitude [$^\circ$E]', fontsize=16)
ax1.set_ylabel('Time [model years]', fontsize=16)
ax1.set_title('a) Temperature anomalies', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
cbar.ax.tick_params(labelsize=13)

CS1 = ax2.contourf(salt_lon_70_58S_anomaly, cmap='BrBG_r', levels=np.linspace(-0.001, 0.001, 21), extend='both')
cbar1 = fig.colorbar(CS1, ax=ax2, ticks=[-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.1])  # Set specific ticks
cbar1.ax.tick_params(labelsize=13)
cbar1.set_label(label = 'Salinity anomaly [g/kg]', fontsize=14)
ax2.set_ylabel('Longitude [$^\circ$E]')
#ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_title('b) Salinity anomalies', fontsize=17)
ax2.tick_params(axis='both', which='major', labelsize=14)
plt.tight_layout()
plt.savefig(directory_figures +'SALT_TEMP_anomalies_lon_70_58S_HR_pop.pdf')
plt.show()
# %%
