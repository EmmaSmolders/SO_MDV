#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Nov  4 13:29:19 2024

@author: 6008399

Hovmoller diagrams of the SOM in HR-POP

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

#%% Read in data 

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

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

plt.figure()
plt.plot(time, som_star)
plt.axhline(y=3, color='r', linestyle='--')
plt.xlim(500, 600)
plt.ylim(3.6, 4.4)
plt.xlabel('Time')
plt.ylabel('SOM*')
plt.title('SOM* vs Time')
plt.show()

#%% AMOC strength plot

fig, axs = plt.subplots(1, 1, figsize=(8, 6))

axs.plot(time*0.0003, amoc_strength, color = 'black', linewidth=1.5)
axs.set_ylabel('Volume transport [Sv]', fontsize = 18)
axs.set_xlabel('Freshwater flux forcing F$_H$ [Sv]', fontsize=18)
axs.set_title('AMOC strength at 26$^\circ$N', fontsize=22)
axs.set_ylim(-1,25)
axs.set_xlim(0,0.18)
axs.set_xticks([0, 0.06, 0.12, 0.18])
axs.tick_params(axis='both', which='major', labelsize=16)
axs.grid()

ax3 	= fig.add_axes([0.05, 0.155, 0.5, 0.4], projection = ccrs.Orthographic(-30, 10))

ax3.coastlines(resolution='110m')
ax3.gridlines()
ax3.add_feature(cfeature.LAND, zorder=10)
ax3.set_global()


lon     = np.arange(0, 361)
lat     = np.arange(-90, 91)
field   = np.ones((len(lat), len(lon))) * -0.35
CS      = ax3.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon     = np.arange(-100, -5)
lat     = np.arange(20, 43)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax3.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

lon     = np.arange(-100, 3)
lat     = np.arange(42, 51)
field   = np.ones((len(lat), len(lon))) * 0.35
CS      = ax3.contourf(lon, lat, field, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

ax3.text(320, 38, '$+F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=15, transform=ccrs.PlateCarree())
ax3.text(340, -10, '$-F_H$', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=15, transform=ccrs.PlateCarree())

x_1	= np.arange(-81, -9.99, 0.1)
y_1	= np.zeros(len(x_1)) + 26.0
y_2	= np.arange(24, 28.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(24, 28.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax3.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax3.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

plt.tight_layout()
plt.savefig(directory_figures +'AMOC_strength_HRPOP.png', dpi=1200)
plt.show()

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

ax4 = fig.add_axes([0.72, 0.195, 0.302, 0.150], projection=ccrs.SouthPolarStereo())
ax4.set_extent([-50, 50, -80, -33], crs=ccrs.PlateCarree())

ax4.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
ax4.add_feature(cfeature.OCEAN, zorder=0, edgecolor='black')

gl = ax4.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

#Define region for som regions
lon_min, lon_max = -50, 0
lat_min, lat_max = -50, -35

lon_min_star, lon_max_star = -30, 20
lat_min_star, lat_max_star = -55, -40

#ACC transect
acc_lon_min, acc_lon_max = -68, -64
acc_lat_min, acc_lat_max = -66.5, -55

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



plt.tight_layout()
plt.savefig(directory_figures +'SOM_AMOC_transport_HR_pop.pdf')
plt.show()

#%%

fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))

ax2.plot(time[0:300], som[0:300], label='SOM', color='black', linewidth=2)
ax2.set_title('Southern Ocean Mode (SOM)', fontsize = 20)#' and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_ylabel('SOM index [$^\circ$C]', color='black', fontsize=16)
ax2.tick_params(labelcolor='black', labelsize=14)
ax2.grid(True)

plt.tight_layout()
plt.savefig(directory_figures +'SOM_HR_pop.png', dpi=1200)
plt.show()

#%%

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 9))

ax1.plot(time, amoc_strength, label='AMOC Strength', color='blue')
ax1.set_title('a) AMOC strength at 26$^\circ$N', fontsize=20)
#ax1.set_xlabel('Time [model years]')
ax1.set_ylabel('Volume transport [Sv]', fontsize=16)
ax1.tick_params(labelcolor='black', labelsize=14)
ax1.grid(True)

ax2.plot(time, som, label='SOM', color='black', linewidth=2)
ax2.set_title('b) Southern Ocean Mode (SOM)', fontsize = 20)#' and depth averaged SOM (0-1000m)')
#ax2.set_xlabel('Time [model years]', fontsize=16)
ax2.set_ylabel('SOM index [$^\circ$C]', color='black', fontsize=16)
ax2.tick_params(labelcolor='black', labelsize=14)
ax2.grid(True)

ax3.plot(time, acc_transport, label='ACC transport', color='green', linewidth=2)
ax3.set_title('c) Drake Passage transport', fontsize = 20)#' and depth averaged SOM (0-1000m)')
ax3.set_xlabel('Time [model years]', fontsize=16)
ax3.set_ylabel('Volume transport [Sv]', color='black', fontsize=16)
ax3.tick_params(labelcolor='black', labelsize=14)
ax3.grid(True)

plt.tight_layout()
plt.savefig(directory_figures +'SOM_AMOC_ACC_HR_pop.png', dpi=1200)
plt.show()

#%%

plt.figure(figsize=(8,4))
plt.plot(time, acc_transport, color='green', linewidth=2)
plt.title('Drake Passage transport', fontsize=20)
plt.xlabel('Time [model years]', fontsize=16)
plt.ylabel('Volume transport [Sv]', fontsize=16)
plt.ylim(94,125)
plt.xlim(0,600)
plt.grid(True)
plt.tick_params(labelcolor='black', labelsize=14)

plt.tight_layout()
plt.savefig(directory_figures +'Drake_HR_pop.png', dpi=1200)
plt.show()

#%% Fourier analysis

import numpy as np

# Example time series
t_begin=450
t_end = 600
t = time[t_begin:t_end]
#y = som_star[t_begin:t_end]
#y = MLD_max[t_begin:t_end]
#y = amoc_strength[t_begin:t_end]
y = acc_transport[t_begin:t_end]

plt.figure()
plt.plot(t,y)

# Fourier transform
Y = np.fft.fft(y)
freqs = np.fft.fftfreq(len(y), d=(t[1] - t[0]))

plt.figure()
plt.plot(freqs,Y)

# Only consider the positive frequencies
pos_mask = freqs > 0
Ypos = Y[pos_mask]
freqs_pos = freqs[pos_mask]

plt.figure()
plt.plot(freqs_pos,np.abs(Ypos)**2)

# Dominant frequency (index of max amplitude)
dom_idx = np.argmax(np.abs(Ypos))
dominant_freq = freqs_pos[dom_idx]
dominant_period = 1 / dominant_freq

print("Dominant frequency:", dominant_freq)
print("Dominant period:", dominant_period)

#%%

#Time vector, frequency and amplitude
                    #total time length [kyr]
dt = 1                          #time spacing [kyr]

#Plotting data
figure, axis = plt.subplots(2, 1, figsize=(15,8))
figure.tight_layout()
figure.subplots_adjust(hspace=0.3)
axis[0].set_title('Sine wave with a frequency of 0.01 kyr$^{-1}$ and an amplitude of 1', size = 16) 
axis[0].plot(t, y)
axis[0].set_xlabel('Time [kyr]', size = 14)
axis[0].set_ylabel('Amplitude', size = 14)

# Frequency domain representation
fourierTransform = np.fft.fft(y)
power = np.abs(fourierTransform)**2
frequency = np.fft.fftfreq(int((t_end-t_begin)/dt), dt)
pos_mask = frequency > 0
powerpos = power[pos_mask]
frequencypos = frequency[pos_mask]

# Frequency domain representation
axis[1].set_title('Power of the frequency components', size = 16)
axis[1].plot(abs(frequencypos), powerpos)
axis[1].set_xlabel('Frequency [yr$^{-1}$]', size = 14)
axis[1].set_ylabel('Spectral power', size = 14)
axis[1].set_xlim([0,0.05])

#%%

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

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
lat_min_wgkp, lat_max_wgkp = -90, -50

#SO30
lon_min_so30, lon_max_so30 = -180, 180
lat_min_so30, lat_max_so30 = -90, -30

fig = plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.SouthPolarStereo())

ax.set_extent([-50, 20, -80, -23], crs=ccrs.PlateCarree())

ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
ax.add_feature(cfeature.OCEAN, zorder=0, edgecolor='black')

gl = ax.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False

x_1	= np.arange(lon_min, lon_max + 0.1, 0.1)
y_1	= np.zeros(len(x_1)) + lat_min
y_2	= np.zeros(len(x_1)) + lat_max

y_3 = np.arange(lat_min, lat_max + 0.1, 0.1)
x_2 = np.zeros(len(y_3)) + lon_min
x_3 = np.zeros(len(y_3)) + lon_max

ax.plot(x_1, y_1, '-k', linewidth = 4.0, transform=ccrs.PlateCarree(), zorder=10)
ax.plot(x_1, y_2, '-k', linewidth = 4.0, transform=ccrs.PlateCarree(), zorder=10)
ax.plot(x_2, y_3, '-k', linewidth = 4.0, transform=ccrs.PlateCarree(), zorder=10)
ax.plot(x_3, y_3, '-k', linewidth = 4.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_star	= np.arange(lon_min_star, lon_max_star + 0.1, 0.1)
y_1_star	= np.zeros(len(x_1)) + lat_min_star
y_2_star	= np.zeros(len(x_1)) + lat_max_star

y_3_star = np.arange(lat_min_star, lat_max_star + 0.1, 0.1)
x_2_star = np.zeros(len(y_3)) + lon_min_star
x_3_star = np.zeros(len(y_3)) + lon_max_star

#ax.plot(x_1_star, y_1_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_1_star, y_2_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_2_star, y_3_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_3_star, y_3_star, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_acc	= np.arange(acc_lon_min, acc_lon_max + 0.1, 0.1)
y_1_acc	= np.zeros(len(x_1_acc)) + acc_lat_min
y_2_acc	= np.zeros(len(x_1_acc)) + acc_lat_max

y_3_acc = np.arange(acc_lat_min, acc_lat_max + 0.1, 0.1)
x_2_acc = np.zeros(len(y_3_acc)) + -66
x_3_acc = np.zeros(len(y_3_acc)) + acc_lon_max

#ax.plot(x_1_acc, y_1_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_1_acc, y_2_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_2_acc, y_3_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_3_acc, y_3_acc, '-g', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp


#ax.plot(x_1_wgkp, y_1_wgkp, '-', color='orange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_1_wgkp, y_2_wgkp, '-', color='orange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_2_wgkp, y_3_wgkp, '-', color='orange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax.plot(x_3_wgkp, y_3_wgkp, '-', color='orange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

plt.tight_layout()
plt.savefig(directory_figures +'region_SOM.png', dpi=1200)

#%%

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(10, 8))

# First subplot
ax1.plot(time, amoc_strength, label='AMOC Strength', color='blue')
ax1.set_title('a) AMOC Transport (0-1000m)')
ax1.set_xlabel('Time [model years]')
ax1.set_ylabel('AMOC Strength [Sv]')
ax1.grid(True)

# Second subplot with two y-axes
ax2.plot(time, som, label='SOM', color='black')
ax2.set_title('b) Southern Ocean Mode (SOM) and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('SOM [$^\circ$C]', color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(time, som_star, label='SOM$^*$', color='red')
ax2_twin.set_ylabel('SOM$^*$[$^\circ$C]', color='red')
ax2_twin.tick_params(axis='y', labelcolor='red')

# Third subplot
ax3.plot(time, acc_transport, color='green')
ax3.set_title('c) Drake Passage transport')
ax3.set_xlabel('Time [model years]')
ax3.set_ylabel('Volume transport [Sv]')
ax3.grid(True)

# Fourth subplot (example placeholder)
ax4.plot(time, som_star, label='SOM$^*$', color='purple')
ax4.set_title('d) SOM$^*$')
ax4.set_xlabel('Time [model years]')
ax4.set_ylabel('SOM$^*$[$^\circ$C]')
ax4.grid(True)

plt.tight_layout()
plt.show()

#%%

fh2 = netcdf.Dataset(directory + 'DSL_Hovmoller.nc','r')

lat = fh2.variables['lat'][:] #Latitude [deg N]
DSL = fh2.variables['DSL'][:] #Dynamic sea level [cm]

fh2.close()

fh3 = netcdf.Dataset(directory + 'TEMP_SALT_Hovmoller_depth_300-700m.nc','r')

temp = fh3.variables['TEMP'][:] #Temperature [degC]
salt = fh3.variables['SALT'][:] #Salinity [g/kg]

fh3.close()

DSL_anomaly     = np.zeros((len(time), len(lat)))
temp_anomaly    = np.zeros((len(time), len(lat)))
salt_anomaly    = np.zeros((len(time), len(lat)))

for lat_i in range(len(lat)):    
    DSL_anomaly[:, lat_i] = DSL[:, lat_i] - np.mean(DSL[0:100,lat_i])
    temp_anomaly[:, lat_i] = temp[:, lat_i] - np.mean(temp[0:100,lat_i])
    salt_anomaly[:, lat_i] = salt[:, lat_i] - np.mean(salt[0:100,lat_i])

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

fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(15, 6))

CS = ax1.contourf(time, lat, temp_anomaly.transpose(), cmap='seismic', levels=np.linspace(-5.1, 5.1, 21), extend='both')
cbar = fig.colorbar(CS, label = 'Temperature anomaly [$^\circ$C]')
cbar.set_label(label = 'Temperature anomaly [$^\circ$C]', fontsize=14)
ax1.set_ylabel('Latitude [$^\circ$N]', fontsize=16)
ax1.set_xlabel('Time [model years]', fontsize=16)
ax1.set_title('a) Temperature anomalies', fontsize=17)
ax1.tick_params(axis='both', which='major', labelsize=14)
cbar.ax.tick_params(labelsize=13)

CS1 = ax2.contourf(time, lat, salt_anomaly.transpose(), cmap='seismic', levels=np.linspace(-0.51, 0.51, 21), extend='both')
cbar1 = fig.colorbar(CS1)
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

# Create a 1x2 grid of subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey=True)

# Subplot 1: Temperature anomaly negative SOM* cycle 1
im1 = axes[0].contourf(lat, depth, np.mean(temp_negative_som_1, axis=0) - np.mean(temp_1, axis=0),
                       levels=np.arange(-0.3, 0.31, 0.05), extend='both', cmap='RdBu_r')
axes[0].set_title('Negative SOM* Cycle 1')
axes[0].set_ylabel('Depth [m]')
axes[0].set_xlabel('Latitude')
axes[0].set_ylim(depth[-1], 0)

# Subplot 2: Temperature anomaly positive SOM* cycle 1
im2 = axes[1].contourf(lat, depth, np.mean(temp_positive_som_1, axis=0) - np.mean(temp_1, axis=0),
                       levels=np.arange(-0.3, 0.31, 0.05), extend='both', cmap='RdBu_r')
axes[1].set_title('Positive SOM* Cycle 1')
axes[1].set_xlabel('Latitude')
axes[1].set_ylim(depth[-1], 0)

# Add a colorbar
fig.colorbar(im2, ax=axes, orientation='vertical', fraction=0.02, pad=0.04, label='Temperature Anomaly')

# Adjust layout
plt.tight_layout()
plt.show()

# %%
