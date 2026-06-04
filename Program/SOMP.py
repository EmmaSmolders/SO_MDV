#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 12:17:03 2025

@author: 6008399

SOM-P

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

#%%

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/HR_POP/Zenodo/Data_final/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%% Read in data 

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

#fh = netcdf.Dataset(directory + 'SOM_P_lat_-60--45_lon_175-250.nc','r')

#time        = fh.variables['time'][:] #model years
#som_p         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

#fh.close()

fh = netcdf.Dataset(directory + 'SOM_P_lat_-60--45_lon_175-210.nc','r')

time        = fh.variables['time'][:] #model years
som_p         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

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

#%%

plt.figure()
plt.plot(time, som, label='SOM-A')
plt.plot(time, som_p, label='SOM-P')
plt.plot(time, som_p_v2, label='SOM-P v2')
plt.legend()
plt.axvline(x=63)
plt.axvline(x=113)
plt.axvline(x=324)
plt.axvline(x=378)
plt.axvline(x=500)
plt.axvline(x=600)
plt.axvline(x=410)
plt.axvline(x=480)

#%% Plot results

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8))

ax1.plot(time, amoc_strength, label='AMOC Strength', color='darkgreen')
ax1.set_title('a) AMOC strength at 26$^\circ$N', fontsize=15)
#ax1.set_xlabel('Time [model years]')
ax1.set_ylabel('Volume transport [Sv]', fontsize=13)
ax1.tick_params(labelcolor='black')
ax1.grid(True)
#ax1.set_xlim(0,0.18)
ax1.set_xlim(0,600)
ax1.tick_params(axis='both', labelsize=11)

ax2.plot(time, som, label='SOM', color='black')
ax2.set_title('b) Southern Ocean Mode (SOM)', fontsize=15)#' and depth averaged SOM (0-1000m)')
#ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('SOM [$^\circ$C]', color='black', fontsize=13)
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)


ax2_twin = ax2.twinx()
#ax2_twin.plot(time*0.0003, som_star, label='SOM$^*$', color='red')
ax2_twin.plot(time, som_p, label='SOM-P', color='blue')
#ax2_twin.set_ylabel('SOM$^*$[$^\circ$C]', color='red')
ax2_twin.set_ylabel('SOM-P [$^\circ$C]', color='blue', fontsize=13)
ax2_twin.tick_params(axis='y', labelcolor='blue')
#ax2.set_xlim(0,0.18)
ax2.set_xlim(0,600)
ax2.tick_params(axis='both', labelsize=11)
ax2_twin.tick_params(axis='both', labelsize=11)

ax3.plot(time, acc_transport, color='red')
ax3.set_title('c) Drake Passage transport', fontsize=15)
#ax3.set_xlabel('Freshwater flux forcing F$_H$ [Sv]')
ax3.set_xlabel('Time [model year]', fontsize=13)
ax3.set_ylabel('Volume transport [Sv]', fontsize=13)
ax3.set_ylim(94,125)
#ax3.set_xlim(0,0.18)
ax3.set_xlim(0,600)
ax3.grid(True)
ax3.tick_params(axis='both', labelsize=11)

ax4 = fig.add_axes([0.7, 0.205, 0.302, 0.150], projection=ccrs.SouthPolarStereo())
ax4.set_extent([-180, 60, -90, -33], crs=ccrs.PlateCarree())

ax4.add_feature(cfeature.LAND, zorder=0, edgecolor='black')
ax4.add_feature(cfeature.OCEAN, zorder=0, edgecolor='black')

gl = ax4.gridlines(draw_labels=False, dms=True, x_inline=False, y_inline=False)
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = False

#Define region for som regions
lon_min, lon_max = -50, 0
lat_min, lat_max = -50, -35

lon_min_star, lon_max_star = 150, 190
lat_min_star, lat_max_star = -70, -60

lon_min_p, lon_max_p = 170, 210
lat_min_p, lat_max_p = -60, -45

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
y_1_star	= np.zeros(len(x_1_star)) + lat_min_star
y_2_star	= np.zeros(len(x_1_star)) + lat_max_star

y_3_star = np.arange(lat_min_star, lat_max_star + 0.1, 0.1)
x_2_star = np.zeros(len(y_3_star)) + lon_min_star
x_3_star = np.zeros(len(y_3_star)) + lon_max_star

#ax4.plot(x_1_star, y_1_star, '-', color='darkorange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax4.plot(x_1_star, y_2_star, '-', color='darkorange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax4.plot(x_2_star, y_3_star, '-', color='darkorange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax4.plot(x_3_star, y_3_star, '-', color='darkorange', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_p	= np.arange(lon_min_p, lon_max_p + 0.1, 0.1)
y_1_p	= np.zeros(len(x_1_p)) + lat_min_p
y_2_p	= np.zeros(len(x_1_p)) + lat_max_p

y_3_p = np.arange(lat_min_p, lat_max_p + 0.1, 0.1)
x_2_p = np.zeros(len(y_3_p)) + lon_min_p
x_3_p = np.zeros(len(y_3_p)) + lon_max_p

ax4.plot(x_1_p, y_1_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_1_p, y_2_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_2_p, y_3_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_3_p, y_3_p, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)


x_1_acc	= np.arange(acc_lon_min, acc_lon_max + 0.1, 0.1)
y_1_acc	= np.zeros(len(x_1_acc)) + acc_lat_min
y_2_acc	= np.zeros(len(x_1_acc)) + acc_lat_max

y_3_acc = np.arange(acc_lat_min, acc_lat_max + 0.1, 0.1)
x_2_acc = np.zeros(len(y_3_acc)) + -66
x_3_acc = np.zeros(len(y_3_acc)) + acc_lon_max

ax4.plot(x_1_acc, y_1_acc, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_1_acc, y_2_acc, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
ax4.plot(x_2_acc, y_3_acc, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

x_1_wgkp	= np.arange(lon_min_wgkp, lon_max_wgkp + 0.1, 0.1)
y_1_wgkp	= np.zeros(len(x_1_wgkp)) + lat_min_wgkp
y_2_wgkp	= np.zeros(len(x_1_wgkp)) + lat_max_wgkp

y_3_wgkp = np.arange(lat_min_wgkp, lat_max_wgkp + 0.1, 0.1)
x_2_wgkp = np.zeros(len(y_3_wgkp)) + lon_min_wgkp
x_3_wgkp = np.zeros(len(y_3_wgkp)) + lon_max_wgkp

#ax4.plot(x_1_wgkp, y_1_wgkp, '-', color='green', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax4.plot(x_1_wgkp, y_2_wgkp, '-', color='green', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax4.plot(x_2_wgkp, y_3_wgkp, '-', color='green', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)
#ax4.plot(x_3_wgkp, y_3_wgkp, '-', color='green', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder=10)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_1_SOM_OS.pdf')
plt.show()

#%% Dominant periods

t = time[450:600]
y = som_p[450:600]

plt.figure()
plt.plot(t, y)

# Fourier transform
Y = np.fft.fft(y)
freqs = np.fft.fftfreq(len(y), d=(t[1] - t[0]))

plt.figure()
plt.plot(freqs, Y)

# Only consider the positive frequencies:
pos_mask = freqs > 0
Ypos = Y[pos_mask]
freqs_pos = freqs[pos_mask]

plt.figure()
plt.plot(freqs_pos, np.abs(Ypos)**2)

# Dominant frequency (index of max amplitude)
sorted_indices = np.argsort(np.abs(Ypos))[::-1]  # Sort indices by amplitude in descending order
dominant_idx = sorted_indices[0]  # Index of the largest amplitude
second_dominant_idx = sorted_indices[1]  # Index of the second largest amplitude

dominant_freq = freqs_pos[dominant_idx]
second_dominant_freq = freqs_pos[second_dominant_idx]

dominant_period = 1 / dominant_freq
second_dominant_period = 1 / second_dominant_freq

print("Dominant frequency:", dominant_freq)
print("Dominant period:", dominant_period)
print("Second dominant frequency:", second_dominant_freq)
print("Second dominant period:", second_dominant_period)

#%% Perform wavelet analysis with statistical significance testing (wavelet is usefull to tell us what periods exist and how they evolve over time)

import pycwt as wavelet

#time series
ts = acc_transport 
dt = 1.0  #1 year

#remove mean and normalize (important!)
ts = (ts - np.mean(ts)) / np.std(ts)

#Estimate AR(1) red noise 
alpha, _, _ = wavelet.ar1(ts)

#omega0 = 6 is the non-dimensional frequency, which controls the balance between time and frequency resolution. A common choice is 6, which provides a good balance for many applications. Higher values give better frequency resolution but worse time resolution, while lower values give better time resolution but worse frequency resolution.
mother = wavelet.Morlet(6)

#Compute CWT
wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(ts, dt, wavelet = mother)

power = np.abs(wave) ** 2
period = 1 / freqs

#Compute local significance levels (is for the 2D wavelet power spectrum, time vs period). Tells us whether the power at each period is significant relative to the red-noise background
signif, fft_theor = wavelet.significance(
    1.0, dt, scales, 0, alpha,
    significance_level=0.95,
    wavelet=mother)

sig95 = signif[:, None]  # expand for broadcasting
significant = power / sig95

global_ws = power.mean(axis=1)

#Compute global significance levels (is for the global wavelet spectrum, i.e. time mean power at each period). Identify dominant significant period overall
dof = ts.size - scales
signif_global, fft_theor = wavelet.significance(
    np.var(ts),
    dt,
    scales,
    sigma_test=1,
    alpha=alpha,
    significance_level=0.95,
    dof=dof,
    wavelet=mother)

peak_idx = np.argmax(global_ws)
print("Dominant period:", period[peak_idx])

sig_idx = np.where(global_ws > signif_global)[0]
print("Significant periods:", period[sig_idx])

time = np.arange(len(ts))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# --- Wavelet power spectrum ---
cf = axes[0].contourf(time, period, power, levels=20, extend='both')
axes[0].contour(time, period, significant, levels=[1], colors='k', linewidths=1, label='95% significance')
axes[0].plot(time, coi, 'w--', linewidth=1, label='Cone of Influence')
axes[0].set_ylim(100, 10)          
axes[0].set_yscale('linear')      
# axes[0].invert_yaxis()           
axes[0].axvline(415, color='white', linewidth=2, label='Onset AMOC collapse')
axes[0].set_xlabel("Time")
axes[0].set_ylabel("Period (years)")
axes[0].set_title("Wavelet Power Spectrum (ACC transport)")

fig.colorbar(cf, ax=axes[0], label="Power")

axes[0].legend(loc='lower right')

# --- Global wavelet spectrum ---
axes[1].plot(global_ws, period, label='Global wavelet spectrum')
axes[1].plot(signif_global, period, '--', label='95% significance')
axes[1].set_yscale('log')
axes[1].invert_yaxis()
axes[1].set_xlabel("Power")
axes[1].set_ylabel("Period")
axes[1].set_title("Global Wavelet Spectrum")
axes[1].legend()

plt.tight_layout()
plt.show()

#%%

def compute_wavelet(ts, dt=1.0, omega0=6):
    # remove mean and normalize
    ts = np.asarray(ts)
    ts = (ts - np.mean(ts)) / np.std(ts)

    # AR1 red-noise estimate
    alpha, _, _ = wavelet.ar1(ts)

    # wavelet
    mother = wavelet.Morlet(omega0)

    # CWT
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(ts, dt, wavelet=mother)

    power = np.abs(wave) ** 2
    period = 1 / freqs

    # local significance
    signif, fft_theor = wavelet.significance(
        1.0, dt, scales, sigma_test=0, alpha=alpha,
        significance_level=0.95,
        wavelet=mother
    )

    sig95 = signif[:, None]
    significant = power / sig95

    return {
        "power": power,
        "period": period,
        "coi": coi,
        "significant": significant,
        "time": np.arange(len(ts))
    }


# --- your 4 time series ---
series_dict = {
    "ACC transport": acc_transport,
    "SOM-P": som_p,
    "SOM": som,
    "AMOC strength": amoc_strength
}

results = {}
for name, ts in series_dict.items():
    results[name] = compute_wavelet(ts, dt=1.0, omega0=6)


fig, axes = plt.subplots(2, 2, figsize=(12, 8))
axes = axes.reshape(2, 2)

# common color levels
all_power = np.concatenate([r["power"].ravel() for r in results.values()])
vmin = np.nanpercentile(all_power, 5)
vmax = np.nanpercentile(all_power, 95)
levels = np.linspace(vmin, vmax, 20)

for i in range(2):
    for j in range(2):
        ax = axes[i, j]
        name = list(results.keys())[i*2 + j]
        r = results[name]

        cf = ax.contourf(
            r["time"], r["period"], r["power"],
            levels=levels, extend="both")

        ax.contour(
            r["time"], r["period"], r["significant"],
            levels=[1], colors="k", linewidths=1)

        ax.plot(r["time"], r["coi"], "w--", linewidth=1, label='Cone of Influence')

        ax.set_ylim(100, 0)
        ax.set_yscale("linear")

        ax.axvline(415, color="white", linewidth=2, label='Onset AMOC collapse')

        # --- labels control ---
        if j == 0:  # left column
            ax.set_ylabel("Period [years]", fontsize=11)
        else:
            ax.set_ylabel("")

        if i == 1:  # bottom row
            ax.set_xlabel("Time [model years]", fontsize=11)
        else:
            ax.set_xlabel("")

        ax.set_title(name)

        # --- individual colorbar ---
        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label("Power")

        ax.legend(loc="upper left")

plt.tight_layout()
#plt.savefig(directory_figures + 'SOM_wavelet.pdf')
plt.show()

#%%

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
axes = np.array(axes)

# select only these three
selected_results = {
    "a) ACC transport": results["ACC transport"],
    "b) SOM": results["SOM"],
    "c) SOM-P": results["SOM-P"]
}

# common color levels
all_power = np.concatenate([r["power"].ravel() for r in selected_results.values()])
vmin = np.nanpercentile(all_power, 5)
vmax = np.nanpercentile(all_power, 95)
levels = np.linspace(0, 30, 21)

for k, (ax, (name, r)) in enumerate(zip(axes, selected_results.items())):
    cf = ax.contourf(
        r["time"], r["period"], r["power"],
        levels=levels, extend="max")

    ax.contour(
        r["time"], r["period"], r["significant"],
        levels=[1], colors="k", linewidths=1)

    ax.plot(r["time"], r["coi"], "w--", linewidth=1, label='Cone of Influence')
    ax.axvline(415, color="darkred", linewidth=2, label='Onset AMOC collapse')

    ax.set_ylim(100, 5)
    ax.set_xlim(0, 601)
    ax.set_yscale("linear")

    if k == 0:
        ax.set_ylabel("Period [years]", fontsize=11)
        ax.legend(loc="upper left")
    else:
        ax.set_ylabel("")

    ax.set_xlabel("Time [model years]", fontsize=11)
    ax.set_title(name, fontsize=12)

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Power")

plt.tight_layout()
plt.savefig(directory_figures + 'SOM_wavelet.pdf')
plt.show()

#%%

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'MLD_year_1-600_Australia.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_AU= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_year_1-600_area_averaged_WGKP.nc','r')

time_temp_wgkp = fh.variables['time'][:] #Time SOM cycle 1
depth_temp_wgkp= fh.variables['depth'][:] 
temp_wgkp = fh.variables['TEMP'][:]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/MLD_year_1-600_NZ.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_NZ = fh.variables['MLD_max'][:] 
MLD = fh.variables['MLD'][:]

fh.close()

plt.figure(figsize=(8,4))
plt.contourf(time_temp_wgkp, depth_temp_wgkp, (temp_wgkp - np.mean(temp_wgkp[0:50], axis=0)).transpose(), levels=np.linspace(-0.4, .4, 21), extend='both', cmap = 'RdBu_r')
plt.colorbar(label='Temperature anomaly [$^\circ$C]')
plt.plot(time_mld, MLD_max, color='black')
#plt.vlines(500, 5000, 0, color='black', linestyle='--')
#plt.vlines(400, 5000, 0, color='black', linestyle='--')
plt.title('MLD maximum and temperature anomaly (WGKP)')
plt.ylim(5000, 0)
plt.ylabel('Depth [m]')
plt.xlabel('Time [model years]')
#plt.savefig(directory_figures +'MLD_WGKP_max_1-600_HR_pop.pdf')

#%%

fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))


ax2.plot(time, som_p, label='SOM-P', color='black')
ax2.set_title('SOM-P and MXL WGKP')#' and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('SOM-P [$^\circ$C]', color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(time_mld, MLD_max, label='MLD', color='red')
ax2_twin.set_ylabel('MLD maximum WGKP [m]]]', color='red')
ax2_twin.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(0,600)

plt.savefig(directory_figures + 'SOMP_max_MLD_WGKP.pdf')

#%%

def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))

ax2.plot(time, som_p, label='SOM-P', color='black')
ax2.set_title('SOM-P and MXL AU')#' and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('SOM-P [$^\circ$C]', color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(Moving_average(MLD_max_AU, 20), label='MLD', color='red')
ax2_twin.set_ylabel('MLD maximum AU [m]]]', color='red')
ax2_twin.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(0,600)

plt.savefig(directory_figures + 'SOMP_max_MLD_AU.pdf')

#%%
fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))


ax2.plot(time, acc_transport, label='ACC transport', color='green')
ax2.set_title('ACC transport and MXL WGKP and NZ')#' and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('ACC transport [Sv]', color='green')
ax2.tick_params(axis='y', labelcolor='green')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(time_mld, MLD_max, label='MLD', color='red')
ax2_twin.plot(time_mld, MLD_max_NZ, label='max', color='black')
ax2_twin.set_ylabel('MLD maximum [m]', color='red')
ax2_twin.tick_params(axis='y', labelcolor='red')
ax2.set_xlim(0,600)

plt.savefig(directory_figures + 'ACC_max_MLD_WGKP_NZ.pdf')

#%%

fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))


ax2.plot(time, acc_transport, label='ACC transport', color='green')
ax2.set_title('ACC transport and SOM-P')#' and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('ACC transport [Sv]', color='green')
ax2.tick_params(axis='y', labelcolor='green')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(time, som_p, label='SOM-P', color='black')
ax2_twin.set_ylabel('SOM-P [C]', color='black')
ax2_twin.tick_params(axis='y', labelcolor='black')
ax2.set_xlim(0,600)

plt.savefig(directory_figures + 'ACC_SOM_P.pdf')

#%%

fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))


ax2.plot(time_mld, MLD_max_NZ, label='NZ', color='green')
ax2.set_title('MLD NZ and SOM-P')#' and depth averaged SOM (0-1000m)')
ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('MLD NZ [m]', color='green')
ax2.tick_params(axis='y', labelcolor='green')
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(Moving_average(som_p, 20), label='SOM-P', color='black')
ax2_twin.set_ylabel('SOM-P [C]', color='black')
ax2_twin.tick_params(axis='y', labelcolor='black')
ax2.set_xlim(0,600)

#plt.savefig(directory_figures + 'ACC_SOM_P.pdf')


#%%

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data

fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

axs[0].plot(norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), color='red', label='max MLD WGKP')
axs[0].plot(norm(Moving_average(som, 20), Moving_average(som, 20)), color='blue', label='SOM index')
#axs[0].plot(norm(Moving_average(acc_transport,20), Moving_average(acc_transport,20)), color='blue', label='Drake Passage')
axs[0].grid()
axs[0].legend()
axs[0].set_ylabel('Normalised scale')
axs[0].set_title('a) MLD WGKP region and SOM index')
axs[0].set_ylim(-3, 3)

axs[1].plot(norm(Moving_average(MLD_max_NZ, 20), Moving_average(MLD_max_NZ, 20)), color='darkorange', label='max MLD NZ')
axs[1].plot(norm(Moving_average(som_p, 20), Moving_average(som_p, 20)), color='green', label='SOM-P index')
#axs[1].plot(norm(Moving_average(acc_transport,20), Moving_average(acc_transport,20)), color='blue', label='Drake Passage')
axs[1].set_title('b) MLD NZ region and SOM-P index')
plt.grid()
plt.legend()
plt.xlabel('Time [model years]')
plt.ylabel('Normalised scale')
plt.savefig(directory_figures + 'ACC_SOM_P_max_MLD_NZ_WGKP.pdf')

#%%

fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

# Subplot 1: MLD WGKP region and SOM index
ax1 = axs[0]
ax2 = ax1.twinx()  # Create a second y-axis for the first subplot

ax1.plot(norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), color='red', label='max MLD WGKP')
ax2.plot(norm(Moving_average(som, 20), Moving_average(som, 20)), color='blue', label='SOM index')

ax1.set_ylabel('Normalised MLD WGKP', color='red')
ax2.set_ylabel('Normalised SOM index', color='blue')
ax1.tick_params(axis='y', labelcolor='red')
ax2.tick_params(axis='y', labelcolor='blue')

ax1.grid()
ax1.set_title('a) MLD WGKP region and SOM index')

# Subplot 2: MLD NZ region and SOM-P index
ax3 = axs[1]
ax4 = ax3.twinx()  # Create a second y-axis for the second subplot

ax3.plot(norm(Moving_average(MLD_max_NZ, 20), Moving_average(MLD_max_NZ, 20)), color='darkorange', label='max MLD NZ')
ax4.plot(norm(Moving_average(som_p, 20), Moving_average(som_p, 20)), color='green', label='SOM-P index')

ax3.set_ylabel('Normalised MLD NZ', color='darkorange')
ax4.set_ylabel('Normalised SOM-P index', color='green')
ax3.tick_params(axis='y', labelcolor='darkorange')
ax4.tick_params(axis='y', labelcolor='green')

ax3.set_title('b) MLD NZ region and SOM-P index')
ax3.grid()

# Shared x-axis
axs[1].set_xlabel('Time [model years]')

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig(directory_figures + 'ACC_SOM_P_max_MLD_NZ_WGKP.pdf')
plt.show()

#%%

plt.figure(figsize=(8,4))
plt.plot(norm(Moving_average(MLD_max_NZ, 20), Moving_average(MLD_max_NZ, 20)), color='darkorange', label='max MLD NZ')
plt.plot(norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), color='red', label='max MLD WGKP')
plt.plot(norm(Moving_average(acc_transport,20), Moving_average(acc_transport,20)), color='royalblue', label='Drake Passage transport')
plt.grid()
plt.legend()
plt.ylim(-4.2, 4.2)
plt.xlim(0,600)
plt.xlabel('Time [model years]')
plt.ylabel('Normalised scale')
plt.title('MLD NZ, MLD WGKP and Drake Passage transport')
plt.savefig(directory_figures + 'MLD_NZ_WGKP_Drake.pdf')

#%% Read in data 

fh = netcdf.Dataset(directory + 'Ocean/TEMP_SALT_DENS_year_63-114_zonal_averaged_lon_170-290_transect_SO_Pacific.nc','r')
#fh = netcdf.Dataset(directory + 'TEMP_SALT_year_1-51_zonal_averaged_55W_5W_transect_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1      = fh.variables['TEMP'][:]   #temperature
salt_1      = fh.variables['SALT'][:]   #salinity
dens_1      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/TEMP_SALT_DENS_year_500-600_zonal_averaged_lon_170-290_transect_SO_Pacific.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_3      = fh.variables['TEMP'][:]   #temperature
salt_3      = fh.variables['SALT'][:]   #salinity
dens_3      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_63-114_zonal_averaged_lon_170-290_transect_SO_Pacific.nc','r')

u_1      = fh.variables['U_VEL'][:]   #zonal velocity
lat_u_1  = fh.variables['lat'][:]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_500-600_zonal_averaged_lon_170-290_transect_SO_Pacific.nc','r')

u_3      = fh.variables['U_VEL'][:]   #zonal velocity
lat_u_3   = fh.variables['lat'][:]

fh.close()

plt.figure()
plt.contourf(lat_u_3, depth, u_3)

plt.figure()
plt.contourf(lat_u_1, depth, u_1)

#%% Take meridional density gradient at every grid point using a various bin size

bin_size = 20  # number of indices to each side

drho_dy1 = ma.masked_all((len(depth), len(lat)))
drho_dy3 = ma.masked_all((len(depth), len(lat)))

# Convert degrees of latitude to meters
deg2m = 111e3

for depth_i in range(len(depth)):
    for lat_i in range(bin_size, len(lat) - bin_size):

        lat_1 = lat_i - bin_size
        lat_2 = lat_i + bin_size

        del_lat_m = np.abs(lat[lat_2] - lat[lat_1]) * deg2m

        del_dens1 = dens_1[depth_i, lat_2] - dens_1[depth_i, lat_1]
        del_dens3 = dens_3[depth_i, lat_2] - dens_3[depth_i, lat_1]

        drho_dy1[depth_i, lat_i] = del_dens1 / del_lat_m
        drho_dy3[depth_i, lat_i] = del_dens3 / del_lat_m



#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=5)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

#plt.suptitle('SOM cycle 3 - 1, Pacific Ocean (170$^\circ$E - 110$^\circ$W)', fontsize=16)

CS = axs[0,0].contourf(lat, depth, temp_3 - temp_1, levels = np.arange(-3, 3.01, 0.2), extend = 'both',  cmap = 'RdBu_r')#, norm = divnorm, cmap = 'RdBu_r')
#axs[0,0].set_xlim(-75,-10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-3, 3.01, 1))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_3 - salt_1, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
#axs[0,1].set_xlim(-75,-10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy3 - drho_dy1, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
#axs[1,0].set_xlim(-75,-10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS4 = axs[1,1].contourf(lat_u_3, depth, u_3 - u_1, levels = np.arange(-0.05, 0.05, 0.001), extend = 'both', cmap = 'RdBu_r')
#axs[1,1].set_xlim(-75,-10)
axs[1,1].set_ylim(depth[-1], 0)
axs[1,1].set_xticklabels(['70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS4, ticks = np.arange(-0.05, .051, 0.05))
cbar.set_label('Zonal velocity difference [m/s]', fontsize = 11)
axs[1,1].set_title('d) Zonal velocity', fontsize=14)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_S2_SOM_OS.pdf')
plt.show()

#%% Read in data 

fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_63-114_zonal_averaged_lon_50-150_transect_SO_Indian.nc','r')
#fh = netcdf.Dataset(directory + 'TEMP_SALT_year_1-51_zonal_averaged_55W_5W_transect_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1      = fh.variables['TEMP'][:]   #temperature
salt_1      = fh.variables['SALT'][:]   #salinity
dens_1      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_lon_50-150_transect_SO_Indian.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_3      = fh.variables['TEMP'][:]   #temperature
salt_3      = fh.variables['SALT'][:]   #salinity
dens_3      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'UVEL_year_63-114_zonal_averaged_lon_50-150_transect_SO_Indian.nc','r')

u_1      = fh.variables['U_VEL'][:]   #zonal velocity

fh.close()

fh = netcdf.Dataset(directory + 'UVEL_year_324-378_zonal_averaged_lon_50-150_transect_SO_Indian.nc','r')

u_2      = fh.variables['U_VEL'][:]   #zonal velocity

fh.close()

fh = netcdf.Dataset(directory + 'UVEL_year_500-600_zonal_averaged_lon_50-150_transect_SO_Indian.nc','r')

u_3      = fh.variables['U_VEL'][:]   #zonal velocity
lat_u     = fh.variables['lat'][:]

fh.close()

#%% Take meridional density gradient at every grid point using a various bin size

bin_size = 20  # number of indices to each side

drho_dy1 = ma.masked_all((len(depth), len(lat)))
drho_dy3 = ma.masked_all((len(depth), len(lat)))

# Convert degrees of latitude to meters
deg2m = 111e3

for depth_i in range(len(depth)):
    for lat_i in range(bin_size, len(lat) - bin_size):

        lat_1 = lat_i - bin_size
        lat_2 = lat_i + bin_size

        del_lat_m = np.abs(lat[lat_2] - lat[lat_1]) * deg2m

        del_dens1 = dens_1[depth_i, lat_2] - dens_1[depth_i, lat_1]
        del_dens3 = dens_3[depth_i, lat_2] - dens_3[depth_i, lat_1]

        drho_dy1[depth_i, lat_i] = del_dens1 / del_lat_m
        drho_dy3[depth_i, lat_i] = del_dens3 / del_lat_m

#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=5)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

plt.suptitle('SOM cycle 3 - 1, Indian Ocean (50$^\circ$E - 150$^\circ$E)', fontsize=16)

CS = axs[0,0].contourf(lat, depth, temp_3 - temp_1, levels = np.arange(-3, 3.01, 0.2), extend = 'both',  cmap = 'RdBu_r')#, norm = divnorm, cmap = 'RdBu_r')
#axs[0,0].set_xlim(-70,-35)
axs[0,0].set_ylim(depth[-1], 0)
#axs[0,0].set_xticklabels(['70', '65', '60', '55', '50', '45', '40', '35'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-3, 3.01, 1))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_3 - salt_1, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
#axs[0,1].set_xlim(-70,-35)
axs[0,1].set_ylim(depth[-1], 0)
#axs[0,1].set_xticklabels(['70', '65', '60', '55', '50', '45', '40', '35'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy3 - drho_dy1, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
#axs[1,0].set_xlim(-70,-35)
axs[1,0].set_ylim(depth[-1], 0)
#axs[1,0].set_xticklabels(['70', '65', '60', '55', '50', '45', '40', '35'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS4 = axs[1,1].contourf(lat_u, depth, u_3 - u_1, levels = np.arange(-0.05, 0.051, 0.01), extend = 'both', cmap = 'RdBu_r')
#axs[1,1].set_xlim(-70,-35)
axs[1,1].set_ylim(depth[-1], 0)
#axs[1,1].set_xticklabels(['70', '65', '60', '55', '50', '45', '40', '35'])
cbar	= colorbar(CS4, ticks = np.arange(-0.1, .11, 0.05))
cbar.set_label('Zonal velocity difference [m/s]', fontsize = 11)
axs[1,1].set_title('d) Zonal velocity', fontsize=14)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.tight_layout()
plt.savefig(directory_figures +'TEMP_SALT_MERGRAD_UVEL_difference_Indian_3_1_HR_pop.pdf')
plt.show()

#%% test uvel pacific

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_500-600_zonal_averaged_lon_170-290_transect_SO_Pacific.nc','r')

u_3      = fh.variables['U_VEL'][:]   #zonal velocity
depth    = fh.variables['depth'][:]
lat_u_3   = fh.variables['lat'][:]

fh.close()

plt.figure()
plt.contourf(lat_u_3, depth, u_3)



