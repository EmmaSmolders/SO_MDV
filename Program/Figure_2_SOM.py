#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 16:04:55 2025

@author: 6008399

Figure 2 SOM paper

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

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_1-100_zonal_averaged_60W_0W_transect_SO.nc','r')

lat_u_1  = fh.variables['lat'][:]
u_1      = fh.variables['U_VEL'][:]   #zonal velocity

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/UVEL_year_500-600_zonal_averaged_60W_0W_transect_SO.nc','r')

u_2      = fh.variables['U_VEL'][:]   #zonal velocity

fh.close()

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
            
            print(lat[lat_1])
            print(lat[lat_2])
            
            #sys.exit()
            
            if lat[lat_2] < lat[lat_1]:
                print('Something is going wrong here')
                sys.exit()
                print(lat[lat_1])
                print(lat[lat_2])
            del_lat_actual[lat_i] = np.abs(lat[lat_2] - lat[lat_1])
            
            del_dens1 = dens_1[depth_i, lat_2] - dens_1[depth_i, lat_1]
            del_dens2 = dens_2[depth_i, lat_2] - dens_2[depth_i, lat_1]
            drho_dy1[depth_i, lat_i] = del_dens1/del_lat_actual[lat_i]
            drho_dy2[depth_i, lat_i] = del_dens2/del_lat_actual[lat_i]
            
            #sys.exit()
            
#%%

bin_size = 20
drho_dy1 = ma.masked_all((len(depth), len(lat)))
drho_dy2 = ma.masked_all((len(depth), len(lat)))

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
        del_dens2 = dens_2[depth_i, lat_2] - dens_2[depth_i, lat_1]

        drho_dy1[depth_i, lat_i] = del_dens1 / del_lat_m
        drho_dy2[depth_i, lat_i] = del_dens2 / del_lat_m




#%%

divnorm = mcolors.TwoSlopeNorm(vmin=-2, vcenter=0, vmax=6)
divnorm_salt = mcolors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=3)

fig, axs = plt.subplots(2, 2, figsize=(12, 6))

CS = axs[0,0].contourf(lat, depth, temp_2 - temp_1, levels = np.arange(-2, 6.01, 0.2), extend = 'both', norm = divnorm, cmap = 'RdBu_r')
axs[0,0].set_xlim(-75,-10)
axs[0,0].set_ylim(depth[-1], 0)
axs[0,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[0,0].set_ylabel('Depth [m]', fontsize=12)
cbar	= colorbar(CS, ticks = np.arange(-2, 6.01, 2))
cbar.set_label('Temperature difference [$^\circ$C]', fontsize = 11)
axs[0,0].set_title('a) Temperature', fontsize=14)

CS2 = axs[0,1].contourf(lat, depth, salt_2 - salt_1, levels = np.arange(-1, 1.01, 0.1), extend = 'both', cmap = 'BrBG_r')
axs[0,1].set_xlim(-75,-10)
axs[0,1].set_ylim(depth[-1], 0)
axs[0,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS2, ticks = np.arange(-1, 1.01, 0.5))
cbar.set_label('Salinity difference [g/kg]', fontsize = 11)
axs[0,1].set_title('b) Salinity', fontsize=14)

CS3 = axs[1,0].contourf(lat, depth, drho_dy2 - drho_dy1, levels = np.linspace(-0.0000005, 0.0000005, 21), extend = 'both', cmap = 'PuOr_r')

#plt.legend()
cbar	= colorbar(CS3, ticks = np.arange(-0.0000005, 0.00000051, 0.0000005))
cbar.set_label(r'$\Delta \rho / \Delta y$ [kg/m$^4$]', fontsize = 11)
axs[1,0].set_title('c) Merdidional density gradient', fontsize=14)
axs[1,0].set_xlim(-75,-10)
axs[1,0].set_ylim(depth[-1], 0)
axs[1,0].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
axs[1,0].set_ylabel('Depth [m]', fontsize=12)
axs[1,0].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

CS4 = axs[1,1].contourf(lat_u_1, depth, u_2 - u_1, levels = np.arange(-0.1, 0.11, 0.01), extend = 'both', cmap = 'RdBu_r')
axs[1,1].set_xlim(-75,-10)
axs[1,1].set_ylim(depth[-1], 0)
axs[1,1].set_xticklabels(['80', '70', '60', '50', '40', '30', '20', '10'])
cbar	= colorbar(CS4, ticks = np.arange(-0.1, .11, 0.05))
cbar.set_label('Zonal velocity difference [m/s]', fontsize = 11)
axs[1,1].set_title('d) Zonal velocity', fontsize=14)
axs[1,1].set_xlabel('Latitude [$^\circ$S]', fontsize=12)

plt.tight_layout()
plt.savefig(directory_figures +'Figure_3_SOM_OS.pdf')
plt.show()