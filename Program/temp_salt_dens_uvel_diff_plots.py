#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 17:45:09 2024

@author: 6008399

Temperature, salinity, density and zonal velocity differences HR pop in SOM region

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
#import gsw

#Making pathway to folder with all data
directory = r'/Users/6008399/Documents/PhD/HR_POP/netcdf/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%% Read in data 

fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_45-91_zonal_averaged_60W_0W_transect_SO.nc','r')
#fh = netcdf.Dataset(directory + 'TEMP_SALT_year_1-51_zonal_averaged_55W_5W_transect_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1      = fh.variables['TEMP'][:]   #temperature
salt_1      = fh.variables['SALT'][:]   #salinity
dens_1      = fh.variables['PD'][:]     #potential density

fh.close()

#fh = netcdf.Dataset(directory + 'DENS_year_1-51_zonal_averaged_55W_5W_transect_SO.nc','r')

#depth       = fh.variables['depth'][:]  #depth
#lat         = fh.variables['lat'][:]    #latitude
#dens_1      = fh.variables['PD'][:]     #potential density

#fh.close()

#fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_550-600_zonal_averaged_55W_5W_transect_SO.nc','r')
fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_495-571_zonal_averaged_60W_0W_transect_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_2      = fh.variables['TEMP'][:]   #temperature
salt_2      = fh.variables['SALT'][:]   #salinity
dens_2      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory + 'UVEL_year_45-91_zonal_averaged_60W_0W_transect_SO.nc','r')

u_1      = fh.variables['U_VEL'][:]   #zonal velocity

fh.close()

fh = netcdf.Dataset(directory + 'UVEL_year_495-571_zonal_averaged_60W_0W_transect_SO.nc','r')

u_2      = fh.variables['U_VEL'][:]   #zonal velocity

fh.close()

#%% Read in data

fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_63-114_zonal_averaged_60W_0W_transect_SO.nc','r')
#fh = netcdf.Dataset(directory + 'TEMP_SALT_year_1-51_zonal_averaged_55W_5W_transect_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
dy          = fh.variables['DYT'][:]    #meridional grid spacing
temp_1      = fh.variables['TEMP'][:]   #temperature
salt_1      = fh.variables['SALT'][:]   #salinity
dens_1      = fh.variables['PD'][:]     #potential density

fh.close()

#fh = netcdf.Dataset(directory + 'DENS_year_1-51_zonal_averaged_55W_5W_transect_SO.nc','r')

#depth       = fh.variables['depth'][:]  #depth
#lat         = fh.variables['lat'][:]    #latitude
#dens_1      = fh.variables['PD'][:]     #potential density

#fh.close()

#fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_550-600_zonal_averaged_55W_5W_transect_SO.nc','r')
fh = netcdf.Dataset(directory + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_60W_0W_transect_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_2      = fh.variables['TEMP'][:]   #temperature
salt_2      = fh.variables['SALT'][:]   #salinity
dens_2      = fh.variables['PD'][:]     #potential density

fh.close()

plt.figure()
plt.contourf(lat, depth, temp_2 - temp_1)
plt.colorbar()

#%% Take meridional density gradient at every grid point using a various bin size

bin_size = 20
drho_dy1 = ma.masked_all((len(depth), len(lat)))
drho_dy2 = ma.masked_all((len(depth), len(lat)))

del_lat = np.abs(lat[bin_size - bin_size] - lat[bin_size + bin_size])
del_lat_actual = ma.masked_all(len(lat))

print(del_lat)

for depth_i in range(len(depth)):
    print(depth_i)
    for lat_i in range(len(lat) - bin_size):
        if lat_i >= bin_size:
            lat_1 = np.abs(lat - (lat[lat_i] - del_lat/2)).argmin()
            lat_2 = np.abs(lat - (lat[lat_i] + del_lat/2)).argmin()
            if lat_1 > lat_2:
                print('Something is going wrong here')
                sys.exit()
                print(lat_1)
                print(lat_2)
            del_lat_actual[lat_i] = np.abs(lat[lat_2] - lat[lat_1])
            
            del_dens1 = dens_1[depth_i, lat_1] - dens_1[depth_i, lat_2]
            del_dens2 = dens_2[depth_i, lat_1] - dens_2[depth_i, lat_2]
            drho_dy1[depth_i, lat_i] = del_dens1/del_lat
            drho_dy2[depth_i, lat_i] = del_dens2/del_lat
            
            #sys.exit()


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


#%%

plt.figure()
plt.contour(lat, depth, dens_1, levels = [1027.1], colors = 'magenta', linewidths = 1, label='0-50')
plt.contour(lat, depth, dens_2, levels = [1027.1], linestyles = '--', colors = 'magenta', linewidths = 1, label='550-600')
plt.contour(lat, depth, dens_1, levels = [1027.3], colors = 'k', linewidths = 1, label='0-50')
plt.contour(lat, depth, dens_2, levels = [1027.3], linestyles = '--', colors = 'k', linewidths = 1, label='550-600')
plt.contour(lat, depth, dens_1, levels = [1027.5], colors = 'g', linewidths = 1, label='0-50')
plt.contour(lat, depth, dens_2, levels = [1027.5], linestyles = '--', colors = 'g', linewidths = 1, label='550-600')
plt.contour(lat, depth, dens_1, levels = [1027.8], colors = 'b', linewidths = 1, label='0-50')
plt.contour(lat, depth, dens_2, levels = [1027.8], linestyles = '--', colors = 'b', linewidths = 1, label='550-600')

plt.ylim(depth[-1], 0)

#%% 

plt.figure()
plt.plot(np.mean(dens_1, axis=1) - np.mean(dens_2, axis=1), depth, label='Before - after')
plt.legend()
plt.ylim(depth[-1], 0)
plt.title('Area-averaged PD (75$^\circ$S-10$^\circ$N, 60-0$^\circ$W)')
plt.savefig(directory_figures +'dens_before_after_SOM_HR_pop_areamean_75S_10N_6W_0W.pdf')

plt.figure()
plt.plot(np.mean(dens_1[:,594:725], axis=1) - np.mean(dens_2[:,594:725], axis=1), depth, label='Before-after')
plt.legend()
plt.ylim(depth[-1], 0)
plt.title('Area-averaged PD (45-35$^\circ$S, 60-0$^\circ$W)')
plt.savefig(directory_figures +'dens_before_after_SOM_HR_pop_areamean_45S_35S_6W_0W.pdf')

#%% N^2

dz = np.gradient(depth)  # spacing between adjacent depths

# centered difference
density_before = np.mean(dens_1, axis=1)
drho_dz_before = np.zeros_like(density_before)
drho_dz_before[1:-1] = (density_before[2:] - density_before[:-2]) / (depth[2:] - depth[:-2])

density_after = np.mean(dens_2, axis=1)
drho_dz_after = np.zeros_like(density_after)
drho_dz_after[1:-1] = (density_after[2:] - density_after[:-2]) / (depth[2:] - depth[:-2])

# forward/backward difference at boundaries (if needed)
drho_dz_before[0]  = (density_before[1] - density_before[0]) / (depth[1] - depth[0])
drho_dz_before[-1] = (density_before[-1] - density_before[-2]) / (depth[-1] - depth[-2])

drho_dz_after[0]  = (density_after[1] - density_after[0]) / (depth[1] - depth[0])
drho_dz_after[-1] = (density_after[-1] - density_after[-2]) / (depth[-1] - depth[-2])

# 3) Compute N² = −(g / ρ) (dρ/dz)
g = 9.81  # m/s^2
N2_before = - (g / density_before) * drho_dz_before
N2_after = - (g / density_after) * drho_dz_after

plt.figure()
plt.plot(-N2_before - -N2_after, depth, label='Before - after')
plt.vlines(x=0, ymin=5000, ymax=0, color='black')
plt.title('Area-averaged N$^2$ (75$^\circ$S-10$^\circ$N, 60-0$^\circ$W)')
plt.legend(loc=4)
plt.ylim(2000, 10)
plt.savefig(directory_figures +'N2_before_after_SOM_HR_pop_areamean_75S_10N_6W_0W.pdf')

#%%

lat1 = 594
lat2 = 725

lat3=0
lat4=100

lat5=300
lat6 = 400

# centered difference
density_before = np.mean(dens_1[:,lat1:lat2], axis=1)
drho_dz_before = np.zeros_like(density_before)
drho_dz_before[1:-1] = (density_before[2:] - density_before[:-2]) / (depth[2:] - depth[:-2])

density_after = np.mean(dens_2[:,lat1:lat2], axis=1)
drho_dz_after = np.zeros_like(density_after)
drho_dz_after[1:-1] = (density_after[2:] - density_after[:-2]) / (depth[2:] - depth[:-2])

# forward/backward difference at boundaries (if needed)
drho_dz_before[0]  = (density_before[1] - density_before[0]) / (depth[1] - depth[0])
drho_dz_before[-1] = (density_before[-1] - density_before[-2]) / (depth[-1] - depth[-2])

drho_dz_after[0]  = (density_after[1] - density_after[0]) / (depth[1] - depth[0])
drho_dz_after[-1] = (density_after[-1] - density_after[-2]) / (depth[-1] - depth[-2])

# 3) Compute N² = −(g / ρ) (dρ/dz)
g = 9.81  # m/s^2
N2_before = - (g / density_before) * drho_dz_before
N2_after = - (g / density_after) * drho_dz_after

plt.figure()
plt.plot(-N2_after - -N2_before, depth, label='After-before')
plt.vlines(x=0, ymin=5000, ymax=0, color='black')
plt.title('Area-averaged N$^2$ (45-35$^\circ$S, 60-0$^\circ$W)')
plt.legend(loc=3)
plt.ylim(5000, 10)
plt.xlim(-1e-5, 5e-6)
#plt.savefig(directory_figures +'N2_before_after_SOM_HR_pop_areamean_45S_35S_60W_0W.pdf')

#%%


# Latitude ranges
lat3, lat4 = 100, 330
lat5, lat6 = 430, 594
lat1, lat2 = 594, 725

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

# Compute N² for each latitude range
N2_before_1, N2_after_1 = compute_N2(np.mean(dens_1[:, lat1:lat2], axis=1), np.mean(dens_2[:, lat1:lat2], axis=1), depth)
N2_before_2, N2_after_2 = compute_N2(np.mean(dens_1[:, lat3:lat4], axis=1), np.mean(dens_2[:, lat3:lat4], axis=1), depth)
N2_before_3, N2_after_3 = compute_N2(np.mean(dens_1[:, lat5:lat6], axis=1), np.mean(dens_2[:, lat5:lat6], axis=1), depth)

plt.figure()
plt.plot(-N2_after_2 - -N2_before_2, depth, label='Latitude '+str(int(lat[lat3]))+' - '+str(int(lat[lat4]))+'$^\circ$N')
plt.plot(-N2_after_3 - -N2_before_3, depth, label='Latitude '+str(int(lat[lat5]))+' - '+str(int(lat[lat6]))+'$^\circ$N')
plt.plot(-N2_after_1 - -N2_before_1, depth, label='Latitude '+str(int(lat[lat1]))+' - '+str(int(lat[lat2]))+'$^\circ$N')
plt.vlines(x=0, ymin=5000, ymax=0, color='black')
plt.title('Area-averaged N$^2$ After minus Before (60-0$^\circ$W)')
plt.legend(loc=3)
plt.ylim(5000, 10)
plt.xlim(-1e-5, 1e-5)
plt.savefig(directory_figures +'N2_before_after_SOM_HR_pop_areamean_60-0W.pdf')

#%%
N2_before_wgkp_1_2, N2_after_wgkp_1_2 = compute_N2(np.mean(dens_1_wgkp, axis=0), np.mean(dens_2_wgkp, axis=0), depth)
N2_before_wgkp_1_3, N2_after_wgkp_1_3 = compute_N2(np.mean(dens_1_wgkp, axis=0), np.mean(dens_3_wgkp, axis=0), depth)

plt.figure()
plt.plot(-N2_after_wgkp_1_2 - -N2_before_wgkp_1_2, depth, label='2 vs 1')
plt.plot(-N2_after_wgkp_1_3 - -N2_before_wgkp_1_3, depth, label='3 vs 1')
plt.vlines(x=0, ymin=5000, ymax=0, color='black', linestyle = '--')
plt.title('Area-averaged N$^2$ difference (WGKP)')
plt.xlabel('N$^2$ difference')
plt.ylabel('Depth [m]')
plt.legend(loc=3)
plt.ylim(5000, 10)
plt.xlim(-1.5e-6, 5e-6)
plt.savefig(directory_figures +'N2_before_after_SOM_HR_pop_areamean_WGKP.pdf')

#%%

SA_before = salt_1
pt_before = temp_1
CT_before = gsw.CT_from_pt(SA_before,pt_before)

SA_after = salt_2
pt_after = temp_2
CT_after = gsw.CT_from_pt(SA_after,pt_after)

N2_before1 = ma.masked_all((len(depth)-1, len(lat)))
p_mid_before = ma.masked_all((len(depth)-1, len(lat)))
N2_after1 = ma.masked_all((len(depth)-1, len(lat)))
p_mid_after = ma.masked_all((len(depth)-1, len(lat)))

for lat_i in range(len(lat)):
    N2_before1[:, lat_i], p_mid_before[:, lat_i] = gsw.Nsquared(SA_before[:,lat_i], CT_before[:,lat_i], gsw.p_from_z(-depth, lat[lat_i]), lat[lat_i])
    N2_after1[:, lat_i], p_mid_after[:, lat_i] = gsw.Nsquared(SA_after[:,lat_i], CT_after[:,lat_i], gsw.p_from_z(-depth, lat[lat_i]), lat[lat_i])

#N2, p_mid = gsw.Nsquared(SA[:,400], CT[:,400], gsw.p_from_z(-depth, lat[400]), lat[400])

plt.figure()
plt.plot(np.mean(N2_before1, axis=1), np.mean(p_mid_before, axis=1))
plt.plot(-N2_before, depth)
plt.ylim(3000, 0)

plt.figure()
plt.plot(np.mean(N2_after1, axis=1), np.mean(p_mid_after, axis=1))
plt.plot(-N2_after, depth)
plt.ylim(3000, 0)

plt.figure()
plt.plot(np.mean(N2_before1, axis=1), np.mean(p_mid_before, axis=1), label='Before')
plt.plot(np.mean(N2_after1, axis=1), np.mean(p_mid_after, axis=1), label='After')
plt.legend(loc=4)
plt.ylim(2000, 0)

#%%
plt.figure()
plt.plot(dens_1[:,200], depth, label='Before')
plt.plot(dens_2[:,200], depth, label='After')
plt.legend()
plt.ylim(depth[-1], 0)
plt.savefig(directory_figures +'dens_before_after_SOM_HR_pop2.pdf')

plt.figure()
plt.plot(dens_1[:,300], depth, label='Before')
plt.plot(dens_2[:,300], depth, label='After')
plt.legend()
plt.ylim(depth[-1], 0)
plt.savefig(directory_figures +'dens_before_after_SOM_HR_pop3.pdf')

#%%

lat_idx = 110 #latitude of 0.05 degN (equator)
depth_idx0 = 18 #depth of 465m
depth_idx1 = 21 #depth of 918m
depth_idx2 = 23 #depth of 1380m
depth_idx3 = 26 #depth of 2125m

plt.figure(figsize=(12,6))
plt.contourf(lat, depth, dens_2 - dens_1, levels = np.arange(-0.5, 0.51, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.colorbar()
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx0, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2, levels = [dens_2[depth_idx0, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx1, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2, levels = [dens_2[depth_idx1, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx2, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2, levels = [dens_2[depth_idx2, lat_idx]], linestyles = '--', colors='black')
plt.contour(lat, depth, dens_1, levels = [dens_1[depth_idx3, lat_idx]], colors = 'black')
plt.contour(lat, depth, dens_2, levels = [dens_2[depth_idx3, lat_idx]], linestyles = '--', colors='black')
plt.ylim(depth[-1],0)

#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=6)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2 - temp_1, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')
axs[0,0].set_xlim(-75,10)
axs[0,0].set_ylim(depth[-1], 0)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
#axs[0,0].contour(lat, depth, temp_1, levels = [20], colors = 'k', linewidths = 1)
#axs[0,0].contour(lat, depth, temp_1, levels = [15], colors = 'r', linewidths = 1)
#axs[0,0].contour(lat, depth, temp_1, levels = [8], colors = 'b', linewidths = 1)
#axs[0,0].contour(lat, depth, temp_1, levels = [3], colors = 'k', linewidths = 1)
cbar.set_label('Temperature diff [$^\circ$C]', fontsize = 9)
axs[0,0].set_title('a) Temperature')

CS2 = axs[0,1].contourf(lat, depth, salt_2 - salt_1, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,10)
axs[0,1].set_ylim(depth[-1], 0)
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity diff [g/kg]', fontsize = 9)
axs[0,1].set_title('b) Salinity')

CS3 = axs[1,0].contourf(lat, depth, dens_2 - dens_1, levels = np.arange(-0.5, 0.51, 0.1), extend = 'both', cmap = 'PuOr_r')

axs[1,0].set_xlim(-75,10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].contour(lat, depth, dens_1, levels = [dens_1[depth_idx0, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2, levels = [dens_2[depth_idx0, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1, levels = [dens_1[depth_idx1, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2, levels = [dens_2[depth_idx1, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1, levels = [dens_1[depth_idx2, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2, levels = [dens_2[depth_idx2, lat_idx]], linestyles = '--', colors='black')
axs[1,0].contour(lat, depth, dens_1, levels = [dens_1[depth_idx3, lat_idx]], colors = 'black')
axs[1,0].contour(lat, depth, dens_2, levels = [dens_2[depth_idx3, lat_idx]], linestyles = '--', colors='black')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.5, .51, 0.5))
cbar.set_label('Density diff [kg/m$^3$]', fontsize = 9)
axs[1,0].set_title('c) Potential density')

#CS4 = axs[1,1].contourf(lat, depth, u_2 - u_1, levels = np.arange(-0.1, 0.11, 0.01), extend = 'both', cmap = 'RdBu_r')
#axs[1,1].set_xlim(-75,10)
#axs[1,1].set_ylim(depth[-1], 0)
#cbar	= colorbar(CS4, ticks = np.arange(-0.1, .11, 0.05))
#cbar.set_label('Zonal velocity diff [m/s]', fontsize = 9)
#axs[1,1].set_title('d) Zonal velocity')

plt.tight_layout()
#plt.savefig(directory_figures +'TEMP_SALT_DENS_UVEL_difference_SOM_HR_pop.pdf')
plt.show()

#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=6)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

plt.suptitle('SOM cycle 3-1, zonal averaged')

CS = axs[0,0].contourf(lat, depth, temp_2 - temp_1, levels = np.linspace(-1, 1.5, 11), extend = 'both', norm = divnorm, cmap = 'seismic')
axs[0,0].set_xlim(-78,-50)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-1, 1.51, 0.5))
#axs[0,0].contour(lat, depth, temp_1, levels = [20], colors = 'k', linewidths = 1)
#axs[0,0].contour(lat, depth, temp_1, levels = [15], colors = 'r', linewidths = 1)
#axs[0,0].contour(lat, depth, temp_1, levels = [8], colors = 'b', linewidths = 1)
#axs[0,0].contour(lat, depth, temp_1, levels = [3], colors = 'k', linewidths = 1)
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_2 - salt_1, levels = np.arange(-0.5, .51, 0.05), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-78,-50)
axs[0,1].set_ylim(depth[-1], 0)
cbar	= colorbar(CS2, ticks = np.arange(-0.5, .51, 0.2))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy2 - drho_dy1, levels = np.arange(-0.05, 0.051, 0.001), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.05, 0.0511, 0.05))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
axs[1,0].set_xlim(-78,-50)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$N]', fontsize=12)

#CS4 = axs[1,1].contourf(lat, depth, u_2 - u_1, levels = np.arange(-0.1, 0.11, 0.01), extend = 'both', cmap = 'RdBu_r')
#axs[1,1].set_xlim(-75,10)
#axs[1,1].set_ylim(depth[-1], 0)
#cbar	= colorbar(CS4, ticks = np.arange(-0.1, .11, 0.05))
#cbar.set_label('Zonal velocity difference [m/s]', fontsize = 11)
#axs[1,1].set_title('d) Zonal velocity', fontsize=14)
#axs[1,1].set_xlabel('Latitude [$^\circ$N]', fontsize=12)

plt.tight_layout()
#plt.savefig(directory_figures +'TEMP_SALT_MERGRAD_UVEL_difference_SOM_HR_pop.pdf')
plt.show()

#%% Drake passage transport

fh = netcdf.Dataset(directory + 'Drake_Passage_transport.nc','r')

time                = fh.variables['time'][:]         #time [model years]
transport           = fh.variables['Transport'][:]    #Transport [Sv]

fh.close()

plt.figure(figsize=(8,4))
plt.plot(time, transport)
plt.ylabel('Volume transport [Sv]')
plt.grid()
plt.xlabel('Time [model years]')
plt.title('Volume transport through Drake Passage')

#%%

fig, axs = plt.subplots(1, 2, figsize=(12, 6))

CS = axs[0].contourf(lat, depth, temp_1, levels = np.arange(-2, 2.01, 0.1), extend = 'both', cmap = 'RdBu_r')
axs[0].set_xlim(-75,-55)
axs[0].set_ylim(depth[-1], 0)
#cbar	= colorbar(CS, ticks = np.arange(-6, 6.01, 2))
#cbar.set_label('Temperature [$^\circ$C]', fontsize = 9)
axs[0].set_title('a) Temperature (1-50)')

CS2 = axs[1].contourf(lat, depth, temp_2, levels = np.arange(-2, 2.01, 0.1), extend = 'both',  cmap = 'RdBu_r')
axs[1].set_xlim(-75,-55)
axs[1].set_ylim(depth[-1], 0)
cbar	= colorbar(CS2, ticks = np.arange(-2, 2.01, 1))
cbar.set_label('Temperature [$^\circ$C]', fontsize = 9)
axs[1].set_title('b) Temperature (551-600)')

plt.tight_layout()
plt.savefig(directory_figures +'TEMP_first_last_50y_SOM_HR_pop.pdf')
plt.show()

#%%

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_2 - temp_1, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-50, -35)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature difference (SOM after -  before collapse)')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_SOM2_minus_SOM1_HR_pop.pdf')

#%% Positive and negative SOM phases

fh = netcdf.Dataset(directory + 'TEMP_positive_SOM_star_0-100_years_zonal_averaged_60W_0W_transect_SO.nc','r')

temp_pos_1 = fh.variables['TEMP'][:]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_positive_SOM_star_500-600_years_zonal_averaged_60W_0W_transect_SO.nc','r')

temp_pos_2 = fh.variables['TEMP'][:]

fh.close()

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_pos_2 - temp_pos_1, levels = np.arange(-1., 1.01, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, 0)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature difference (positive SOM* after -  before)')
plt.colorbar()
plt.savefig(directory_figures +'TEMP_diff_posSOMstar2_minus_posSOMstar1_HR_pop.pdf')

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_pos_1 - temp_1, levels = np.arange(-0.5, .51, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, 0)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly positive SOM* before')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_posSOMstar2_minus_posSOMstar1_HR_pop.pdf')

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_pos_2 - temp_2, levels = np.arange(-0.5, .51, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, 0)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly positive SOM* after')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_posSOMstar2_minus_posSOMstar1_HR_pop.pdf')

#%%
fh = netcdf.Dataset(directory + 'TEMP_negative_SOM_star_0-100_years_zonal_averaged_60W_0W_transect_SO.nc','r')

temp_neg_1 = fh.variables['TEMP'][:]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_negative_SOM_star_500-600_years_zonal_averaged_60W_0W_transect_SO.nc','r')

temp_neg_2 = fh.variables['TEMP'][:]

fh.close()

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_neg_2 - temp_neg_1, levels = np.arange(-1., 1.01, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature difference (negative SOM* after -  before)')
plt.colorbar()
plt.savefig(directory_figures +'TEMP_diff_negSOMstar2_minus_negSOMstar1_HR_pop.pdf')

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_neg_1 - temp_1, levels = np.arange(-0.5, .51, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, 0)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly negative SOM* before')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_posSOMstar2_minus_posSOMstar1_HR_pop.pdf')

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, temp_neg_2 - temp_2, levels = np.arange(-0.5, .51, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, 0)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly negative SOM* after')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_posSOMstar2_minus_posSOMstar1_HR_pop.pdf')

#%% Determine positive and negative SOM phases

fh = netcdf.Dataset(directory + 'TEMP_63-114_years_zonal_averaged_60W_0W_transect_SO.nc','r')

time_1 = fh.variables['time'][:]
lat = fh.variables['lat'][:]
depth = fh.variables['depth'][:]
temp_1 = fh.variables['TEMP'][:]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_324-378_years_zonal_averaged_60W_0W_transect_SO.nc','r')

time_2 = fh.variables['time'][:]
temp_2 = fh.variables['TEMP'][:]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_500-600_years_zonal_averaged_60W_0W_transect_SO.nc','r')

time_3 = fh.variables['time'][:]
temp_3 = fh.variables['TEMP'][:]

fh.close()

som_star_1 = som_star[62:114]
som_star_2 = som_star[323:378]
som_star_3 = som_star[499:600]

#Determine positve and negative SOM 1 period
mean_som_1 = np.mean(som_star_1)
std_som_1 = np.std(som_star_1)

#Select the years where the data is higher than the threshold
positive_years_1 = time_1[som_star_1 > mean_som_1 + std_som_1]
positive_som_values_1 = som_star_1[som_star_1 > mean_som_1 + std_som_1]
temp_positive_som_1 = temp_1[som_star_1 > mean_som_1 + std_som_1]

#Select the years where the data is lower than the threshold
negative_years_1 = time_1[som_star_1 < mean_som_1 - std_som_1]
negative_som_values_1 = som_star_1[som_star_1 < mean_som_1 - std_som_1]
temp_negative_som_1 = temp_1[som_star_1 < mean_som_1 - std_som_1]

plt.figure()
plt.title('SOM cycle 1')
plt.plot(time_1, som_star_1)
plt.axvline(x=positive_years_1[0], color='red')
plt.axvline(x=positive_years_1[6], color='red')
plt.axvline(x=positive_years_1[7], color='red')
plt.axvline(x=positive_years_1[-1], color='red')
plt.axvline(x=negative_years_1[0], color='blue')
plt.axvline(x=negative_years_1[-1], color='blue')

#Determine positve and negative SOM 2 period
mean_som_2 = np.mean(som_star_2)
std_som_2 = np.std(som_star_2)

#Select the years where the data is higher than the threshold
positive_years_2 = time_2[som_star_2 > mean_som_2 + std_som_2]
positive_som_values_2 = som_star_2[som_star_2 > mean_som_2 + std_som_2]
temp_positive_som_2 = temp_2[som_star_2 > mean_som_2 + std_som_2]

#Select the years where the data is lower than the threshold
negative_years_2 = time_2[som_star_2 < mean_som_2 - std_som_2]
negative_som_values_2 = som_star_2[som_star_2 < mean_som_2 - std_som_2]
temp_negative_som_2 = temp_2[som_star_2 < mean_som_2 - std_som_2]

plt.figure()
plt.title('SOM cycle 2')
plt.plot(time_2, som_star_2)
plt.axvline(x=positive_years_2[0], color='red')
plt.axvline(x=positive_years_2[-1], color='red')
plt.axvline(x=negative_years_2[0], color='blue')
plt.axvline(x=negative_years_2[-1], color='blue')

#Determine positve and negative SOM 3 period
mean_som_3 = np.mean(som_star_3)
std_som_3 = np.std(som_star_3)

#Select the years where the data is higher than the threshold
#positive_years_3 = time_3[som_star_3 > mean_som_3 + std_som_3]
positive_years_3 = time_3[[15,16, 17, 18, 19,20,21,22,23,24,25,26,27,28,29,30,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]]
positive_som_values_3 = som_star_3[[15,16, 17, 18, 19,20,21,22,23,24,25,26,27,28,29,30,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]]
temp_positive_som_3 = temp_3[[15,16, 17, 18, 19,20,21,22,23,24,25,26,27,28,29,30,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]]

#Select the years where the data is lower than the threshold
negative_years_3 = time_3[som_star_3 < mean_som_3 - std_som_3]
negative_som_values_3 = som_star_3[som_star_3 < mean_som_3 - std_som_3]
temp_negative_som_3 = temp_3[som_star_3 < mean_som_3 - std_som_3]

plt.figure()
plt.title('SOM cycle 3')
plt.plot(time_3, som_star_3)
plt.axhline(y=mean_som_3)
plt.axvline(x=515, color='red')
plt.axvline(x=530, color='red')
plt.axvline(x=586, color='red')
plt.axvline(x=600, color='red')
plt.axvline(x=negative_years_3[0], color='blue')
plt.axvline(x=negative_years_3[2], color='blue')
plt.axvline(x=negative_years_3[3], color='blue')
plt.axvline(x=negative_years_3[-1], color='blue')

#%%

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_negative_som_1, axis=0) - np.mean(temp_1, axis=0), levels = np.arange(-0.3, .31, 0.05), extend = 'both', cmap = 'RdBu_r')
#plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly negative SOM* cycle 1')
plt.colorbar()

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_positive_som_1, axis=0) - np.mean(temp_1, axis=0), levels = np.arange(-.3, .31, 0.05), extend = 'both', cmap = 'RdBu_r')
#plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly positive SOM* cycle 1')
plt.colorbar()

#%%

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_negative_som_2, axis=0) - np.mean(temp_2, axis=0), levels = np.arange(-0.3, .31, 0.05), extend = 'both', cmap = 'RdBu_r')
#plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly negative SOM* cycle 2')
plt.colorbar()

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_positive_som_2, axis=0) - np.mean(temp_2, axis=0), levels = np.arange(-.3, .31, 0.05), extend = 'both', cmap = 'RdBu_r')
#plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly positive SOM* cycle 2')
plt.colorbar()

#%%

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_negative_som_3, axis=0) - np.mean(temp_3, axis=0), levels = np.arange(-.3, .31, 0.05), extend = 'both', cmap = 'RdBu_r')
#plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly negative SOM* cycle 3')
plt.colorbar()

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_positive_som_3, axis=0) - np.mean(temp_3, axis=0), levels = np.arange(-.3, .31, 0.05), extend = 'both', cmap = 'RdBu_r')
#plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature anomaly positive SOM* cycle 3')
plt.colorbar()

#%%

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_negative_som_3, axis=0) - np.mean(temp_negative_som_1, axis=0), levels = np.arange(-1., 1.01, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature difference (negative SOM* after -  before)')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_negSOMstar3_minus_negSOMstar1_HR_pop.pdf')

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_positive_som_3, axis=0) - np.mean(temp_positive_som_1, axis=0), levels = np.arange(-1., 1.01, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature difference (positive SOM* after -  before)')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_posSOMstar3_minus_posSOMstar1_HR_pop.pdf')

plt.figure(figsize = (8,4))
plt.contourf(lat, depth, np.mean(temp_3, axis=0) - np.mean(temp_1, axis=0), levels = np.arange(-1., 1.01, 0.1), extend = 'both', cmap = 'RdBu_r')
plt.xlim(-75, -30)
plt.ylim(depth[-1], 0)
plt.ylabel('Depth [m]')
plt.title('Temperature difference (SOM* after -  before)')
plt.colorbar()
#plt.savefig(directory_figures +'TEMP_diff_SOMstar3_minus_SOMstar1_HR_pop.pdf')


