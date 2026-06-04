#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: 6008399

Wavelet analysis for SOM and SOM-related quantities 

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
directory = r'/Users/6008399/Documents/PhD/HR_POP/netcdf/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%%

def TrendRemover(time, data, trend_type):
	"""Removes trend of choice"""
	
	rank = polyfit(time, data, trend_type)
	fitting = 0.0 
		
	for rank_i in range(len(rank)):
			
		fitting += rank[rank_i] * (time**(len(rank) - 1 - rank_i))

	data -= fitting
	
	return data

#%% Read in data 

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_P_lat_-60--45_lon_175-250.nc','r')

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

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'MLD_year_1-600_Australia.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_AU= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/MLD_year_1-600_NZ.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max_NZ = fh.variables['MLD_max'][:] 
MLD = fh.variables['MLD'][:]

fh.close()

fh = netcdf.Dataset(directory + 'MLD_year_1-600_SO.nc','r')

time_mld        = fh.variables['time'][:]       #time [model years]
MLD_SO_PA       = fh.variables['MLD'][:, 200:600,0:500]  
lat         = fh.variables['lat'][:]    #latitude
lon         = fh.variables['lon'][:]    #latitude
lon_PA          = lon[0:500]
lat_PA          = lat[200:600]

fh.close()


#%% Perform wavelet analysis with statistical significance testing (wavelet is usefull to tell us what periods exist and how they evolve over time)

import pycwt as wavelet

def Moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

#time series
t_start = 450
t_end = 600
ts = TrendRemover(time[t_start:t_end].copy(), amoc_strength[t_start:t_end].copy(), 2)
#ts = Moving_average(ts, n=3)
dt = 1.0  #1 year   

plt.figure()
plt.plot(time[t_start:t_end], ts)

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

#time = np.arange(len(ts))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# --- Wavelet power spectrum ---
cf = axes[0].contourf(time[t_start:t_end], period, power, levels=np.linspace(0,10,20), extend='max')
axes[0].contour(time[t_start:t_end], period, significant, levels=[1], colors='k', linewidths=1, label='95% significance')
axes[0].plot(time[t_start:t_end], coi, 'w--', linewidth=1, label='Cone of Influence')
axes[0].set_ylim(100, 10)          
axes[0].set_yscale('linear')      
# axes[0].invert_yaxis()           
#axes[0].axvline(415, color='white', linewidth=2, label='Onset AMOC collap
# se')
axes[0].set_xlabel("Time [model years]")
axes[0].set_ylabel("Period [years]")
axes[0].set_title("b) Wavelet Power Spectrum (AMOC strength, detrend 2)")

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
plt.savefig(directory_figures + 'AMOC_wavelet_time_'+str(t_start)+'_'+str(t_end)+'.pdf')
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

x = Moving_average(som_p, 4)

plt.figure()
plt.plot(x)

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

#%% MLD wavelet spectra

series_dict = {
    "a) Maximum MLD NZ": MLD_max_NZ,
    "b) Maximum MLD WGKP": MLD_max,
    "c) Maximum MLD AU": MLD_max_AU,
    "d) Maximum MLD PA": np.max(MLD_SO_PA, axis=(1, 2))
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

        ax.axvline(415, color="darkred", linewidth=2, label='Onset AMOC collapse')

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
plt.savefig(directory_figures + 'MLD_wavelet.pdf')
plt.show()
# %%

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data

plt.figure()
plt.plot(norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), label='MLD max WGKP')
plt.plot(norm(Moving_average(MLD_max_AU, 20), Moving_average(MLD_max_AU, 20)), label='MLD max AU')
plt.plot(norm(Moving_average(MLD_max_NZ, 20), Moving_average(MLD_max_NZ, 20)), label='MLD max NZ')
plt.plot(norm(Moving_average(np.max(MLD_SO_PA, axis=(1, 2)), 20), Moving_average(np.max(MLD_SO_PA, axis=(1, 2)), 20)), label='MLD max PA')
#plt.axvline(415, color='darkred', linewidth=2, label='Onset AMOC collapse')
plt.xlabel('Time [model years]')
plt.ylabel('Normalised')
plt.title('Maximum MLD in different regions')
plt.legend()
plt.grid()
plt.savefig(directory_figures + 'MLD_max_time_series.pdf')
plt.show()

# %%

plt.figure()
plt.plot(norm(Moving_average(MLD_max, 20), Moving_average(MLD_max, 20)), label='MLD max WGKP')
plt.plot(norm(Moving_average(acc_transport, 20), Moving_average(acc_transport, 20)), label='ACC transport')
plt.plot(norm(Moving_average(MLD_max_NZ, 20), Moving_average(MLD_max_NZ, 20)), label='MLD max NZ')
#plt.plot(norm(Moving_average(np.max(MLD_SO_PA, axis=(1, 2)), 20), Moving_average(np.max(MLD_SO_PA, axis=(1, 2)), 20)), label='MLD max PA')
#plt.axvline(415, color='darkred', linewidth=2, label='Onset AMOC collapse')
plt.xlabel('Time [model years]')
plt.ylabel('Normalised')
plt.title('Maximum MLD in and ACC transport')
plt.legend()
plt.savefig(directory_figures + 'MLD_NZ_WGKP_ACC_time_series.pdf')
plt.show()
# %%
