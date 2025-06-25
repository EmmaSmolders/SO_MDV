#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 15:22:53 2025

@author: 6008399

Figure S2 SOM paper

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

#%%

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_depth_0-1000m.nc','r')

time               = fh.variables['time'][:]         #time [model years]
som_star           = fh.variables['SOM'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_63-114_window_5_volume_integrated_SO30_finalcode.nc','r')
#fh = netcdf.Dataset(directory + 'TKE_year_63-114_window_5_volume_integrated_SO30.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_1_SO30 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_1_SO30 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_1_SO30 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_324-378_window_5_volume_integrated_SO30_finalcode.nc','r')
#fh = netcdf.Dataset(directory + 'TKE_year_324-378_window_5_volume_integrated_SO30.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 2
TKE_2_SO30 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_2_SO30 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_2_SO30 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_500-600_window_5_volume_integrated_SO30_finalcode.nc','r')
#fh = netcdf.Dataset(directory + 'TKE_year_500-600_window_5_volume_integrated_SO30.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 3
TKE_3_SO30 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_3_SO30 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_3_SO30 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

#fh = netcdf.Dataset(directory + 'KE_year_500-600_window_5_volume_integrated_SO30_fastcode.nc','r')
fh = netcdf.Dataset(directory + 'TKE_year_500-600_window_5_volume_integrated_SO30.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 3
TKE_3_SO30_Rene = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'APE_year_63-114_window_5_volume_integrated_SO30.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_1_SO30 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]
#MPE_1_SO30 = fh.variables['MPE'][:] #Volume integrated mean kinetic energy [J]
#EPE_1_SO30 = fh.variables['EPE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'APE_year_324-378_window_5_volume_integrated_SO30.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_2_SO30 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]
#MPE_2_SO30 = fh.variables['MPE'][:] #Volume integrated mean kinetic energy [J]
#EPE_2_SO30 = fh.variables['EPE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'APE_year_500-600_window_5_volume_integrated_SO30.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
APE_3_SO30 = fh.variables['APE'][:] #Volume integrated total kinetic energy [J]
#MPE_2_SO30 = fh.variables['MPE'][:] #Volume integrated mean kinetic energy [J]
#EPE_2_SO30 = fh.variables['EPE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'Generation_KE_year_63-114_window_5_area_integrated_SO30.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
gKe_1_SO30 = fh.variables['gKe'][:] 
gKm_1_SO30 = fh.variables['gKm'][:]
gKt_1_SO30 = fh.variables['gKtotal'][:]

fh.close()


fh = netcdf.Dataset(directory + 'Generation_KE_year_324-378_window_5_area_integrated_SO30.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 2
gKe_2_SO30 = fh.variables['gKe'][:] 
gKm_2_SO30 = fh.variables['gKm'][:]
gKt_2_SO30 = fh.variables['gKtotal'][:]

fh.close()

fh = netcdf.Dataset(directory + 'Generation_KE_year_500-600_window_5_area_integrated_SO30.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 3
gKe_3_SO30 = fh.variables['gKe'][:] 
gKm_3_SO30 = fh.variables['gKm'][:]
gKt_3_SO30 = fh.variables['gKtotal'][:]

fh.close()

fh = netcdf.Dataset(directory + 'Conversion_PE_KE_year_63-114_5_volume_integrated_SO30_testjuling_Rene.nc','r')

time_1 = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_total_1_SO30 = fh.variables['C_total'][:] #Volume integrated total kinetic energy [J]
C_mean_1_SO30 = fh.variables['C_mean'][:] #Volume integrated mean kinetic energy [J]
C_eddy_1_SO30 = fh.variables['C_eddy'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Conversion_PE_KE_year_324-378_5_volume_integrated_SO30_testjuling_Rene.nc','r')

time_2 = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_total_2_SO30 = fh.variables['C_total'][:] #Volume integrated total kinetic energy [J]
C_mean_2_SO30 = fh.variables['C_mean'][:] #Volume integrated mean kinetic energy [J]
C_eddy_2_SO30 = fh.variables['C_eddy'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Conversion_PE_KE_year_500-600_5_volume_integrated_SO30_testjuling_Rene.nc','r')

time_3 = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_total_3_SO30 = fh.variables['C_total'][:] #Volume integrated total kinetic energy [J]
C_mean_3_SO30 = fh.variables['C_mean'][:] #Volume integrated mean kinetic energy [J]
C_eddy_3_SO30 = fh.variables['C_eddy'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

#fh = netcdf.Dataset(directory + 'Nonzonality_year_63-114_volume_integrated.nc','r')
fh = netcdf.Dataset(directory + 'Nonzonality_year_63-114_volume_integrated_SO30AREA_window_5.nc','r')

time_1 = fh.variables['time_int'][:] #Starting year of window SOM cycle 1
V2U2_1_SO30 = fh.variables['V2_U2_int'][:] #Volume integrated total kinetic energy [J]

fh.close()

#fh = netcdf.Dataset(directory + 'Nonzonality_year_324-378_volume_integrated.nc','r')
fh = netcdf.Dataset(directory + 'Nonzonality_year_324-378_volume_integrated_SO30AREA_window_5.nc','r')

time_2 = fh.variables['time_int'][:] #Starting year of window SOM cycle 1
V2U2_2_SO30 = fh.variables['V2_U2_int'][:] #Volume integrated total kinetic energy [J]

fh.close()

#fh = netcdf.Dataset(directory + 'Nonzonality_year_500-600_volume_integrated.nc','r')
fh = netcdf.Dataset(directory + 'Nonzonality_year_500-600_volume_integrated_SO30AREA_window_5.nc','r')

time_3 = fh.variables['time_int'][:] #Starting year of window SOM cycle 1
V2U2_3_SO30 = fh.variables['V2_U2_int'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_63-114_years_volume_averaged_SO30.nc','r')

time_1_temp = fh.variables['time'][:] #Starting year of window SOM cycle 1
temp_1_SO30 = fh.variables['TEMP'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_324-378_years_volume_averaged_SO30.nc','r')

time_2_temp = fh.variables['time'][:] #Starting year of window SOM cycle 1
temp_2_SO30 = fh.variables['TEMP'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'TEMP_500-600_years_volume_averaged_SO30.nc','r')

time_3_temp = fh.variables['time'][:] #Starting year of window SOM cycle 1
temp_3_SO30 = fh.variables['TEMP'][:] #Volume integrated total kinetic energy [J]

fh.close()

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data

#%%

# Create a 3x3 grid of subplots
fig, axes = plt.subplots(2, 3, figsize=(12, 8))

# Flatten the axes array for easier indexing
axes = axes.flatten()

# Plot 1: Non-zonality
axes[0].plot(V2U2_1_SO30, color='mediumblue', label='SOM cycle 1')
axes[0].plot(V2U2_2_SO30, color='darkorange', label='SOM cycle 2')
axes[0].plot(V2U2_3_SO30, color='crimson', label='SOM cycle 3')
axes[0].set_title('a) $\zeta$', fontsize=14)
axes[0].set_ylabel('$\zeta$ [-]', fontsize = 12)
axes[0].legend()

axes[1].plot(np.gradient(APE_1_SO30, time_1 * 365.25 * 24 * 60 * 60), color='mediumblue')
axes[1].plot(np.gradient(APE_2_SO30, time_2 * 365.25 * 24 * 60 * 60), color='darkorange')
axes[1].plot(np.gradient(APE_3_SO30, time_3 * 365.25 * 24 * 60 * 60), color='crimson')
axes[1].set_title('b) d$P$/dt', fontsize=14)
axes[1].set_ylabel('d$P$/dt [J/s]', fontsize = 12)
#axes[1].legend()

# Plot 5: Mean wind input
axes[2].plot(gKm_1_SO30, color='mediumblue', label='SOM cycle 1')
axes[2].plot(gKm_2_SO30, color='darkorange', label='SOM cycle 2')
axes[2].plot(gKm_3_SO30, color='crimson', label='SOM cycle 3')
axes[2].set_title('c) G($K_m$)', fontsize=14)
axes[2].set_ylabel('G($K_m$) [W]', fontsize=12)
#axes[2].legend()

# Plot 4: EKE
axes[3].plot(EKE_1_SO30, color='mediumblue', label='SOM cycle 1')
axes[3].plot(EKE_2_SO30, color='darkorange', label='SOM cycle 2')
axes[3].plot(EKE_3_SO30, color='crimson', label='SOM cycle 3')
axes[3].set_xlabel('Time [model year]', fontsize = 12)
axes[3].set_ylabel('$K_e$ [J]', fontsize = 12)
axes[3].set_title('d) $K_e$', fontsize=14)
#axes[3].legend()

# Plot 4: EKE
axes[4].plot(TKE_1_SO30, color='mediumblue', label='SOM cycle 1')
axes[4].plot(TKE_2_SO30, color='darkorange', label='SOM cycle 2')
axes[4].plot(TKE_3_SO30, color='crimson', label='SOM cycle 3')
axes[4].set_xlabel('Time [model year]', fontsize = 12)
axes[4].set_title('e) $K$', fontsize=14)
axes[4].set_ylabel('$K$ [J]', fontsize = 12)
#axes[4].legend()

# Plot 8: Eddy conversion
axes[5].plot(C_eddy_1_SO30, color='mediumblue', label='SOM cycle 1')
axes[5].plot(C_eddy_2_SO30, color='darkorange', label='SOM cycle 2')
axes[5].plot(C_eddy_3_SO30, color='crimson', label='SOM cycle 3')
axes[5].set_xlabel('Time [model year]', fontsize = 12)
axes[5].set_title('f) C($P_e$,$K_e$)', fontsize=14)
axes[5].set_ylabel('C($P_e$,$K_e$) [W]', fontsize = 12)
#axes[5].legend()

# Adjust layout
plt.tight_layout()
plt.savefig(directory_figures + 'Magnitude_SOM_energetics_SO30.pdf')
plt.show()
