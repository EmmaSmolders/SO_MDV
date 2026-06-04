#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  8 15:00:32 2025

@author: 6008399

Energetics SO30 with final code and from 1 - 600

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

#%% Read in data SO30

fh = netcdf.Dataset(directory + 'SOM.nc','r')

time        = fh.variables['time'][:] #model years
som         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fh = netcdf.Dataset(directory + 'SOM_depth_0-1000m.nc','r')

time               = fh.variables['time'][:]         #time [model years]
som_star           = fh.variables['SOM'][:]    #Transport [Sv]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/SOM_P_lat_-60--45_lon_175-250.nc','r')

time        = fh.variables['time'][:] #model years
som_p         = fh.variables['SOM'][:]  #Southern Ocean Mode (SOM) [degC]

fh.close()

fig, (ax2) = plt.subplots(1, 1, figsize=(8, 4))

ax2.plot(time, som, label='SOM', color='black')
ax2.set_title('b) Southern Ocean Mode (SOM)')#' and depth averaged SOM (0-1000m)')
#ax2.set_xlabel('Time [model years]')
ax2.set_ylabel('SOM [$^\circ$C]', color='black')
ax2.tick_params(axis='y', labelcolor='black')
ax2.grid(True)

ax2_twin = ax2.twinx()
#ax2_twin.plot(time*0.0003, som_star, label='SOM$^*$', color='red')
ax2_twin.plot(time, som_p, label='SOM-P', color='red')
#ax2_twin.set_ylabel('SOM$^*$[$^\circ$C]', color='red')
ax2_twin.set_ylabel('SOM-P [$^\circ$C]', color='red')
ax2_twin.tick_params(axis='y', labelcolor='red')
#ax2.set_xlim(0,0.18)
ax2.set_xlim(0,600)

#%%

fh = netcdf.Dataset(directory + 'Ocean/KE_year_63-600_window_5_volume_integrated_SO30_finalcode.nc','r')

time_63_600 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_SO30_63_600 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_SO30_63_600 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_SO30_63_600 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_1-70_window_5_volume_integrated_SO30_finalcode.nc','r')

time_1_70 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_SO30_1_70 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
EKE_SO30_1_70 = fh.variables['EKE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/KE_year_390-410_window_5_volume_integrated_SO30_finalcode.nc','r')

time_390_410 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_SO30_390_410 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
EKE_SO30_390_410 = fh.variables['EKE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'KE_year_350-400_window_5_volume_integrated_SO30_finalcode.nc','r')

time_350_400 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_SO30_350_400 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
EKE_SO30_350_400 = fh.variables['EKE'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/KE_year_400-600_window_5_volume_integrated_SO30_finalcode.nc','r')

time_400 = fh.variables['time'][:] #Starting year of window SOM cycle 1
TKE_SO30_400 = fh.variables['TKE'][:] #Volume integrated total kinetic energy [J]
MKE_SO30_400 = fh.variables['MKE'][:] #Volume integrated mean kinetic energy [J]
EKE_SO30_400 = fh.variables['EKE'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

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

fh = netcdf.Dataset(directory + 'MLD_max_year_1-600_WGKP.nc','r')

time_mld = fh.variables['time'][:] #Time SOM cycle 1
MLD_max= fh.variables['MLD_max'][:] 

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Generation_KE_year_400-600_window_5_area_integrated_SO30.nc','r')

time_400 = fh.variables['time'][:] #Starting year of window SOM cycle 1
gKe_SO30 = fh.variables['gKe'][:] 
gKm_SO30 = fh.variables['gKm'][:]
gKtotal_SO30 = fh.variables['gKtotal'][:]

fh.close()

fh = netcdf.Dataset(directory + 'Generation_KE_year_1-200_window_5_area_integrated_SO30.nc','r')

time_1_200 = fh.variables['time'][:] #Starting year of window SOM cycle 1
gKm_SO30_1_200 = fh.variables['gKm'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Generation_KE_year_199-401_window_5_area_integrated_SO30.nc','r')

time_199_401 = fh.variables['time'][:] #Starting year of window SOM cycle 1
gKm_SO30_199_401 = fh.variables['gKm'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Generation_KE_year_190-210_window_5_area_integrated_SO30.nc','r')

time_190_210 = fh.variables['time'][:] #Starting year of window SOM cycle 1
gKm_SO30_190_210 = fh.variables['gKm'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Generation_KE_year_380-410_window_5_area_integrated_SO30.nc','r')

time_380_410 = fh.variables['time'][:] #Starting year of window SOM cycle 1
gKm_SO30_380_410 = fh.variables['gKm'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/Nonzonality_year_63-600_volume_integrated_SO30AREA_window_5.nc','r')

time = fh.variables['time_int'][:] #Starting year of window SOM cycle 1
V2U2_63_600 = fh.variables['V2_U2_int'][:] #Volume integrated total kinetic energy [J]
V2U2_63_600 = V2U2_63_600/3.6473335373774424e+16

fh.close()

fh = netcdf.Dataset(directory + 'Nonzonality_year_1-70_volume_integrated_SO30AREA_window_5.nc','r')

V2U2_1_70 = fh.variables['V2_U2_int'][:] #Volume integrated total kinetic energy [J]

#Forgot to divide by np.sum(volume) in the code
V2U2_1_70 = V2U2_1_70/3.6473335373774424e+16

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Nonzonality_year_1-600_volume_integrated_SO30AREA_window_5.nc','r')

time = fh.variables['time_int'][:] #Starting year of window SOM cycle 1
V2U2_1_600 = fh.variables['V2_U2_int'][:] #Volume integrated total kinetic energy [J]

fh.close()

#%%

fh = netcdf.Dataset(directory + 'Ocean/Conversion_PE_KE_year_1-200_5_volume_integrated_SO30_testjuling.nc','r')

time_1_200_ceddy = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_eddy_SO30_1_200 = fh.variables['C_eddy'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Conversion_PE_KE_year_195-400_5_volume_integrated_SO30_testjuling.nc','r')

time_195_400 = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_eddy_SO30_195_400 = fh.variables['C_eddy'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Conversion_PE_KE_year_395-600_5_volume_integrated_SO30_testjuling.nc','r')

time_395_600 = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_eddy_SO30_395_600 = fh.variables['C_eddy'][:] #Volume integrated total kinetic energy [J]

fh.close()

fh = netcdf.Dataset(directory + 'Ocean/Conversion_PE_KE_year_63-600_5_volume_integrated_SO30_testjuling.nc','r')

time = fh.variables['time'][:] #Starting year of window SOM cycle 1
C_total_SO30_old = fh.variables['C_total'][:] #Volume integrated total kinetic energy [J]
C_mean_SO30 = fh.variables['C_mean'][:] #Volume integrated mean kinetic energy [J]
C_eddy_SO30 = fh.variables['C_eddy'][:] #Volume integrated eddy kinetic energy [J]

fh.close()

def norm(data, data_ref):
    
    #Data has mean of 0 and standard deviation of 1
    norm_data = (data - np.mean(data)) / np.std(data_ref)
    
    return norm_data

fh = netcdf.Dataset(directory + 'Drake_Passage_transport.nc','r')

time_transport                = fh.variables['time'][:]         #time [model years]
acc_transport           = fh.variables['Transport'][:]    #Transport [Sv]

fh.close()

#%% Total energetics

fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, figsize=(8, 10))

ax1.plot(time_1_600, APE_SO30/1e18, label='$\\rho_{ref}$ 0-50', color='royalblue')
ax1.plot(time_350_600, APE_SO30_350_600/1e18, color='royalblue')
ax1.plot(time_som_1, APE_SO30_som_1/1e18, label='$\\rho_{ref}$ SOM', color='cyan')
ax1.plot(time_som_2, APE_SO30_som_2/1e18, color='cyan')
ax1.plot(time_som_3, APE_SO30_som_3/1e18, color='cyan')
ax1.set_title('a) Available potential energy (P)', fontsize = 15)
ax1.set_ylabel('P [EJ]', fontsize = 13)
ax1.tick_params(labelcolor='black')
ax1.grid(True)
ax1.set_xlim(0,600)
ax1.set_ylim(125, 190)
ax1.tick_params(labelsize=12)
ax1.legend(loc = 3, fontsize=11)

ax2.plot(time_1_70, TKE_SO30_1_70/1e18, label='K', color='blue')
ax2.plot(time_63_600, TKE_SO30_63_600/1e18, label='K', color='blue')
ax2.plot(time_350_400, TKE_SO30_350_400/1e18, label='K', color='blue')
ax2.plot(time_390_410, TKE_SO30_390_410/1e18, label='K', color='blue')
ax2.plot(time_400, TKE_SO30_400/1e18, label='K', color='blue')
ax2.set_title('b) Total kinetic energy (K)', fontsize = 15)#' and depth averaged SOM (0-1000m)')
ax2.set_ylabel('K [EJ]', color='black', fontsize = 13)
ax2.tick_params(axis='y', labelcolor='black')
ax2.tick_params(labelsize=12)
ax2.grid(True)

ax2_twin = ax2.twinx()
ax2_twin.plot(time_1_70, EKE_SO30_1_70/1e18, label='K', color='grey')
ax2_twin.plot(time_63_600, EKE_SO30_63_600/1e18, label='K$_e$', color='grey')
ax2_twin.plot(time_350_400, EKE_SO30_350_400/1e18, label='K', color='grey')
ax2_twin.plot(time_390_410, EKE_SO30_390_410/1e18, label='K', color='grey')
ax2_twin.plot(time_400, EKE_SO30_400/1e18, label='K$_e$', color='grey')
ax2_twin.set_ylabel('K$_e$ [EJ]', color='grey', fontsize = 13)
ax2_twin.tick_params(axis='y', labelcolor='grey')
ax2_twin.tick_params(labelsize=12)
ax2.set_xlim(0,600)

ax3.plot(time_400, gKm_SO30/1e12, color='red')
ax3.plot(time_1_200, gKm_SO30_1_200/1e12, color='red')
ax3.plot(time_199_401, gKm_SO30_199_401/1e12, color='red')
ax3.plot(time_190_210, gKm_SO30_190_210/1e12, color='red')
ax3.plot(time_380_410, gKm_SO30_380_410/1e12, color='red')
ax3.set_title('c) Mean wind energy input (G(K$_m$))', fontsize = 15)
ax3.set_ylabel('G(K$_m$) [TW]', fontsize = 13)
ax3.set_xlim(0,600)
ax3.tick_params(labelsize=12)
ax3.set_ylim(1.225, 1.3)
ax3.set_yticks([1.25, 1.3]) 
ax3.set_yticklabels(['1.25', '1.30'])  
ax3.grid(True)

#ax4.plot(time, C_eddy_SO30/1e12, label='old', color='purple')
ax4.plot(time_1_200_ceddy, C_eddy_SO30_1_200/1e12, color='magenta')
ax4.plot(time_195_400, C_eddy_SO30_195_400/1e12, color='magenta')
ax4.plot(time_395_600, C_eddy_SO30_395_600/1e12, color='magenta')
ax4.set_title('d) Conversion eddy potential to eddy kinetic energy (C(P$_e$, K$_e$))', fontsize = 15)
ax4.set_ylabel('C(P$_e$, K$_e$) [TW]', fontsize = 13)
ax4.set_xlim(0,600)
ax4.tick_params(labelsize=12)
ax4.grid(True)
#ax4.legend()

#ax5.plot(time, V2U2_63_600 , color='green')
#ax5.plot(time_1_70, V2U2_1_70, color='green')
ax5.plot(time_1_600, norm(V2U2_1_600, V2U2_1_600), color='green')
ax5.set_title('e) Normalised non-zonality parameter ($\zeta$)', fontsize = 15)
ax5.set_xlabel('Time [model years]', fontsize = 13)
ax5.set_ylabel('$\zeta$ [-]', fontsize = 13)
ax5.set_xlim(0,600)
ax5.tick_params(labelsize=12)
ax5.grid(True)

#plt.suptitle('Energetics SO30', fontsize=16)
plt.tight_layout()
plt.savefig(directory_figures +'Figure_4_SOM_OS.pdf')
plt.show()

#%%

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 6))
#plt.plot(time_som, norm(som))
ax1.plot(time_som_1, norm(TKE_SO30_63_600[0:len(time_som_1)], TKE_SO30_63_600[0:len(time_som_1)]), label='K', color='blue', linewidth=3) #time + 2 because of centered window moving average (and time_1 displays starting year of the window)
ax1.plot(time_som_1, norm(APE_SO30_som_1, APE_SO30_som_1), label='P', color='cyan', linewidth=3)
ax1.plot(time_som_1, norm(gKm_SO30_1_200[62:62+len(time_som_1)], gKm_SO30_1_200[62:62+len(time_som_1)]), label='G(K$_m$)', color='red', linewidth=3)
ax1.plot(time_som_1, norm(C_eddy_SO30_1_200[62:62+len(time_som_1)], C_eddy_SO30_1_200[62:62+len(time_som_1)]), label='C(P$_e$, K$_e$)', color='magenta', linewidth=3)
ax1.plot(time_som_1, norm(V2U2_63_600[0:len(time_som_1)], V2U2_63_600[0:len(time_som_1)]), label='$\zeta$', color='green', linewidth=3)
#ax1.plot(time_som_1, norm(som[63:111], som[63:111]), label = 'SOM', color='black', linewidth=3)
ax1.set_ylim(-2,2.2)
ax1.legend(loc=3, fontsize=11)
ax1.set_ylabel('Normalised scale', fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_xlabel('Time [model year]', fontsize=18)
ax1.set_title('a) SOM cycle 1', fontsize=20)
#ax1.grid()

ax1.axvline(x = time_som_1[0] + np.argmin(APE_SO30_som_1), color='grey')
ax1.axvline(x = time_som_1[0] + np.argmin(TKE_SO30_63_600[0:len(time_som_1)]), color='grey')
ax1.axvline(x = time_som_1[0] + np.argmax(APE_SO30_som_1), color='grey')
ax1.axvline(x = time_som_1[0] + np.argmax(TKE_SO30_63_600[0:len(time_som_1)]), color='grey')

ax1.text(68, 1.9, 'C', fontsize=18, color='black', ha='center', va='center')
ax1.text(80, 1.9, 'D', fontsize=18, color='black', ha='center', va='center')
ax1.text(89, 1.9, 'A', fontsize=18, color='black', ha='center', va='center')
ax1.text(100, 1.9, 'B', fontsize=18, color='black', ha='center', va='center')

time_2 = time[410-63:480-63]

ax2.plot(time_som_2, norm(TKE_SO30_400[10:10+len(time_som_2)], TKE_SO30_63_600[0:len(time_som_1)]), label='TKE', color='blue', linewidth=3)
ax2.plot(time_som_2, norm(APE_SO30_som_2, APE_SO30_som_1), label='APE', color='cyan', linewidth=3)
ax2.plot(time_som_2, norm(gKm_SO30[10:10+len(time_som_2)], gKm_SO30_1_200[62:62+len(time_som_1)]), label='Wind', color='red', linewidth=3)
ax2.plot(time_som_2, norm(C_eddy_SO30_395_600[15:15+len(time_som_2)], C_eddy_SO30_1_200[62:62+len(time_som_1)]), label='PE -> KE', color='magenta', linewidth=3)
ax2.plot(time_som_2, norm(V2U2_63_600[410 - 63:410-63 + len(time_som_2)], V2U2_63_600[0:len(time_som_1)]), label='$\zeta$', color='green', linewidth=3)
#ax2.legend()
ax2.set_ylim(-2,2.2)
ax2.set_xlabel('Time [model year]', fontsize=18)
ax2.tick_params(axis='both', which='major', labelsize=16)
#ax2.grid()
ax2.set_title('b) SOM cycle 2', fontsize=20)

ax2.text(416, 1.9, 'C', fontsize=18, color='black', ha='center', va='center')
ax2.text(431, 1.9, 'D', fontsize=18, color='black', ha='center', va='center')
ax2.text(443, 1.9, 'A', fontsize=18, color='black', ha='center', va='center')
ax2.text(454, 1.9, 'B', fontsize=18, color='black', ha='center', va='center')
ax2.text(468, 1.9, 'C', fontsize=18, color='black', ha='center', va='center')


ax2.axvline(x = time_som_2[0] + np.argmin(APE_SO30_som_2), color='grey')
ax2.axvline(x = time_som_2[0] + np.argmin(TKE_SO30_400[10:10+len(time_som_2)]), color='grey')
ax2.axvline(x = 461, color='grey')
ax2.axvline(x = time_som_2[0] + 8 + np.argmax(TKE_SO30_400[10:10+len(time_som_2)-10]), color='grey')

ax3.plot(time_som_3, norm(TKE_SO30_400[100::], TKE_SO30_63_600[0:len(time_som_1)]), color='blue', linewidth=3)
ax3.plot(time_som_3, norm(APE_SO30_som_3, APE_SO30_som_1), color='cyan', linewidth=3)
ax3.plot(time_som_3, norm(gKm_SO30[100::], gKm_SO30_1_200[62:62+len(time_som_1)]),  color='red', linewidth=3)
#ax3.plot(time_som_3, norm(C_eddy_SO30[500-63::], C_eddy_SO30[0:52]), color='purple', linewidth=3) #OLD ONE!!
ax3.plot(time_som_3, norm(C_eddy_SO30_395_600[-100 + 3::], C_eddy_SO30_1_200[62:62+len(time_som_1)]), color='magenta', linewidth=3) 
ax3.plot(time_som_3, norm(V2U2_63_600[500-63::], V2U2_63_600[0:len(time_som_1)]), color='green', linewidth=3)
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
plt.savefig(directory_figures +'Figure_5_SOM_OS.pdf')
plt.show()


#%%

plt.figure(figsize=(8, 6))
plt.plot(time_1_70, norm(TKE_SO30_1_70, TKE_SO30_1_70), label='K', color='blue')
plt.plot(time_63_600, norm(TKE_SO30_63_600, TKE_SO30_1_70), label='K', color='blue')
plt.plot(time_350_400, norm(TKE_SO30_350_400, TKE_SO30_1_70), label='K', color='blue')
plt.plot(time_390_410, norm(TKE_SO30_390_410, TKE_SO30_1_70), label='K', color='blue')
plt.plot(time_400, norm(TKE_SO30_400, TKE_SO30_1_70), label='K', color='blue')
plt.plot(time_1_200_ceddy, norm(C_eddy_SO30_1_200, C_eddy_SO30_1_200), color='magenta')
plt.plot(time_195_400, norm(C_eddy_SO30_195_400, C_eddy_SO30_1_200), color='magenta')
plt.plot(time_395_600, norm(C_eddy_SO30_395_600, C_eddy_SO30_1_200), color='magenta')
plt.title('Total kinetic energy (K) and conversion eddy potential to eddy kinetic energy (C(P$_e$, K$_e$))', fontsize = 15)
plt.ylabel('K [EJ] and C(P$_e$, K$_e$) [TW]', fontsize = 13)
plt.xlabel('Time [model years]', fontsize = 13)
plt.xlim(0,600)
plt.tick_params(labelsize=12)
#plt.xlim(60,90)
plt.grid(True)
plt.tight_layout()







# %%
