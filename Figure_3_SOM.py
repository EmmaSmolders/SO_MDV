#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 16:34:35 2025

@author: 6008399

Figure 3 SOM paper

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

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
#plt.plot(time_som, norm(som))
ax1.plot(time_1+2, norm(TKE_1_SO30, TKE_1_SO30), label='K', color='blue', linewidth=3) #time + 2 because of centered window moving average (and time_1 displays starting year of the window)
ax1.plot(time_1+2, norm(APE_1_SO30, APE_1_SO30), label='P', color='cyan', linewidth=3)
#ax1.plot(time_1_temp, norm(temp_1_SO30, temp_1_SO30), label='Temp', color='grey', linewidth=3)
#ax1.plot(time[63:115], norm(MLD_max[63:115], MLD_max[63:115]), label='max MLD', color='orange')
#ax1.plot(time[63:115], norm(som[63:115], som[63:115]), label = 'SOM', color='black')
ax1.plot(time[63:115], norm(som_star[63:115], som_star[63:115]), label = 'SOM*', color='black', linewidth=3)
ax1.plot(time_1+2, norm(gKm_1_SO30, gKm_1_SO30), label='G(K$_m$)', color='red', linewidth=3)
ax1.plot(time_1+2, norm(C_eddy_1_SO30, C_eddy_1_SO30), label='C(P$_e$, K$_e$)', color='magenta', linewidth=3)
ax1.plot(time_1+2, norm(V2U2_1_SO30, V2U2_1_SO30), label='$\zeta$', color='green', linewidth=3)
#ax1.plot(time_1+2, norm(V2U2_int_1), '--', label='V$^2$/U$^2$', color='black',)
ax1.set_ylim(-2,2.2)
ax1.legend(loc=3, fontsize=11)
ax1.set_ylabel('Normalized scale', fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_xlabel('Time [model year]', fontsize=18)
ax1.set_title('a) SOM cycle 1', fontsize=20)
#ax1.grid()

ax1.axvline(x = time_1[0] + 2 + np.argmin(APE_1_SO30), color='grey')
ax1.axvline(x = time_1[0] + 2 + np.argmin(TKE_1_SO30)+4, color='grey')
ax1.axvline(x = time_1[0] + 2 + np.argmax(APE_1_SO30), color='grey')
ax1.axvline(x = time_1[0] + 2 + np.argmax(TKE_1_SO30), color='grey')

ax1.text(70, 1.9, 'C', fontsize=18, color='black', ha='center', va='center')
ax1.text(84, 1.9, 'D', fontsize=18, color='black', ha='center', va='center')
ax1.text(94, 1.9, 'A', fontsize=18, color='black', ha='center', va='center')
ax1.text(104, 1.9, 'B', fontsize=18, color='black', ha='center', va='center')

ax2.plot(time_2+2, norm(TKE_2_SO30, TKE_1_SO30), label='TKE', color='blue', linewidth=3)
ax2.plot(time_2+2, norm(APE_2_SO30, APE_1_SO30), label='APE', color='cyan', linewidth=3)
#ax2.plot(time_2_temp, norm(temp_2_SO30, temp_1_SO30), label='Temp', color='grey', linewidth=3)
#ax2.plot(time[324:379], norm(MLD_max[324:379], MLD_max[63:115]), label='max MLD', color='orange')
ax2.plot(time[324:379], norm(som_star[324:379], som_star[63:115]), label = 'SOM*', color='black', linewidth=3)
ax2.plot(time_2+2, norm(gKm_2_SO30, gKm_1_SO30), label='Wind', color='red', linewidth=3)
ax2.plot(time_2+2, norm(C_eddy_2_SO30, C_eddy_1_SO30), label='PE -> KE', color='magenta', linewidth=3)
ax2.plot(time_2+2, norm(V2U2_2_SO30, V2U2_1_SO30), label='$\zeta$', color='green', linewidth=3)
#ax2.legend()
ax2.set_ylim(-2,2.2)
ax2.set_xlabel('Time [model year]', fontsize=18)
ax2.tick_params(axis='both', which='major', labelsize=16)
#ax2.grid()
ax2.set_title('b) SOM cycle 2', fontsize=20)

ax2.text(330, 1.9, 'C', fontsize=18, color='black', ha='center', va='center')
ax2.text(343, 1.9, 'D', fontsize=18, color='black', ha='center', va='center')
ax2.text(352, 1.9, 'A', fontsize=18, color='black', ha='center', va='center')
ax2.text(362, 1.9, 'B', fontsize=18, color='black', ha='center', va='center')

ax2.axvline(x = time_2[0] + 2 + np.argmin(APE_2_SO30), color='grey')
ax2.axvline(x = time_2[0] + 2 + np.argmin(TKE_2_SO30), color='grey')
ax2.axvline(x = time_2[0] + 2 + np.argmax(APE_2_SO30), color='grey')
ax2.axvline(x = time_2[0] + 2 + np.argmax(TKE_2_SO30), color='grey')

ax3.plot(time_3+2, norm(TKE_3_SO30, TKE_1_SO30), color='blue', linewidth=3)
ax3.plot(time_3+2, norm(APE_3_SO30, APE_1_SO30), color='cyan', linewidth=3)
#ax3.plot(time_3_temp, norm(temp_3_SO30, temp_1_SO30), label='Temp', color='grey', linewidth=3)
#ax3.plot(time[500:601], norm(MLD_max[500:601], MLD_max[63:115]), label='max MLD', color='orange')
ax3.plot(time[500:601], norm(som_star[500:601], som_star[63:115]), label = 'SOM*', color='black', linewidth=3)
ax3.plot(time_3+2, norm(gKm_3_SO30, gKm_1_SO30),  color='red', linewidth=3)
ax3.plot(time_3+2, norm(C_eddy_3_SO30, C_eddy_1_SO30), color='magenta', linewidth=3)
ax3.plot(time_3+2, norm(V2U2_3_SO30, V2U2_1_SO30), color='green', linewidth=3)
#ax3.legend()
ax3.set_ylim(-2,2.2)
ax3.set_xlabel('Time [model year]', fontsize=18)
ax3.tick_params(axis='both', which='major', labelsize=16)
#ax3.grid()
ax3.set_title('c) SOM cycle 3', fontsize=20)

#ax2.axvline(x = time_3[0] + 2 + np.argmin(APE_1_SO30), color='grey')
#ax3.axvline(x = 571, color='grey')
#ax2.axvline(x = time_3[0] + 2 + np.argmax(APE_1_SO30), color='grey')
#ax3.axvline(x = time_3[0] + 2 + np.argmax(TKE_3_SO30), color='grey')

plt.tight_layout()
plt.savefig(directory_figures +'Energetics_SOM_SO30.pdf')
plt.show()