#Stratification changes in the four convective regions

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
directory_data = r'/Users/6008399/Documents/PhD/HR_POP/Zenodo/Data_Final/'
directory_figures = r'/Users/6008399/Documents/PhD/HR_POP/Figures/'

#%% Functions

def compute_N2_from_profile(density, depth, rho0=1027, g=9.81):
    """
    Compute vertical density gradient (drho/dz) and buoyancy frequency N²
    from a 1D density profile and depth array.

    Compute N² exactly following the user's original code,
    assuming depth increases downward.
    """

    AU = len(depth)
    n0 = np.zeros(AU)

    for i in range(AU):
        if i == 0:
            # surface: (rho0 - rho1)/(z1 - z0)
            n0[i] = (density[i] - density[i+1]) / (depth[i+1] - depth[i])

        elif i == AU - 1:
            # bottom: (rho_{n-2} - rho_{n-1})/(z_{n-1}-z_{n-2})
            n0[i] = (density[i-1] - density[i]) / (depth[i] - depth[i-1])

        else:
            # centered: (rho[i-1] - rho[i+1]) / (z[i+1] - z[i-1])
            n0[i] = (density[i-1] - density[i+1]) / (depth[i+1] - depth[i-1])

    # Your formula: N² = -(g / rho0) * n0
    N2 = - (g / rho0) * n0
    
    if N2.any() < 0:
        print('N2 is unstable!')

    return n0, N2

#%%

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_zonal_averaged_80E-150E_AU_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1_au      = fh.variables['TEMP'][:]   #temperature
salt_1_au      = fh.variables['SALT'][:]   #salinity
dens_1_au      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_80E-150E_AU_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_500_au      = fh.variables['TEMP'][:]   #temperature
salt_500_au      = fh.variables['SALT'][:]   #salinity
dens_500_au      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_zonal_averaged_-110E--60E_PA_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1_pa      = fh.variables['TEMP'][:]   #temperature
salt_1_pa      = fh.variables['SALT'][:]   #salinity
dens_1_pa      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_-110E--60E_PA_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_500_pa      = fh.variables['TEMP'][:]   #temperature
salt_500_pa      = fh.variables['SALT'][:]   #salinity
dens_500_pa      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_zonal_averaged_160E-190E_NZ_-90N-10N_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat_NZ         = fh.variables['lat'][:]    #latitude
temp_1_nz      = fh.variables['TEMP'][:]   #temperature
salt_1_nz      = fh.variables['SALT'][:]   #salinity
dens_1_nz      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_160E-190E_NZ_-90N-10N_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat_NZ         = fh.variables['lat'][:]    #latitude
temp_500_nz      = fh.variables['TEMP'][:]   #temperature
salt_500_nz      = fh.variables['SALT'][:]   #salinity
dens_500_nz      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_zonal_averaged_-35E-80E_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_1_wgkp      = fh.variables['TEMP'][:]   #temperature
salt_1_wgkp      = fh.variables['SALT'][:]   #salinity
dens_1_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_zonal_averaged_-35E-80E_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lat         = fh.variables['lat'][:]    #latitude
temp_500_wgkp      = fh.variables['TEMP'][:]   #temperature
salt_500_wgkp      = fh.variables['SALT'][:]   #salinity
dens_500_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()
# %% Compute N2 for each profile

n0_lat_1_wgkp = ma.masked_all((len(depth), len(lat)))
n0_lat_500_wgkp = ma.masked_all((len(depth), len(lat)))
n0_lat_1_nz = ma.masked_all((len(depth), len(lat_NZ)))
n0_lat_500_nz = ma.masked_all((len(depth), len(lat_NZ)))
n0_lat_1_au = ma.masked_all((len(depth), len(lat)))
n0_lat_500_au = ma.masked_all((len(depth), len(lat)))
n0_lat_1_pa = ma.masked_all((len(depth), len(lat)))
n0_lat_500_pa = ma.masked_all((len(depth), len(lat)))

N2_lat_1_wgkp = ma.masked_all((len(depth), len(lat)))
N2_lat_500_wgkp = ma.masked_all((len(depth), len(lat)))
N2_lat_1_nz = ma.masked_all((len(depth), len(lat_NZ)))
N2_lat_500_nz = ma.masked_all((len(depth), len(lat_NZ)))
N2_lat_1_au = ma.masked_all((len(depth), len(lat)))
N2_lat_500_au = ma.masked_all((len(depth), len(lat)))
N2_lat_1_pa = ma.masked_all((len(depth), len(lat)))
N2_lat_500_pa = ma.masked_all((len(depth), len(lat)))

for lat_i in range(len(lat)):
    n0_lat_1_wgkp[:, lat_i], N2_lat_1_wgkp[:,lat_i] = compute_N2_from_profile(dens_1_wgkp[:,lat_i], depth)
    n0_lat_500_wgkp[:, lat_i], N2_lat_500_wgkp[:,lat_i] = compute_N2_from_profile(dens_500_wgkp[:,lat_i], depth)
    n0_lat_1_nz[:, lat_i], N2_lat_1_nz[:,lat_i] = compute_N2_from_profile(dens_1_nz[:,lat_i], depth)
    n0_lat_500_nz[:, lat_i], N2_lat_500_nz[:,lat_i] = compute_N2_from_profile(dens_500_nz[:,lat_i], depth)
    n0_lat_1_au[:, lat_i], N2_lat_1_au[:,lat_i] = compute_N2_from_profile(dens_1_au[:,lat_i], depth)
    n0_lat_500_au[:, lat_i], N2_lat_500_au[:,lat_i] = compute_N2_from_profile(dens_500_au[:,lat_i], depth)
    n0_lat_1_pa[:, lat_i], N2_lat_1_pa[:,lat_i] = compute_N2_from_profile(dens_1_pa[:,lat_i], depth)
    n0_lat_500_pa[:, lat_i], N2_lat_500_pa[:,lat_i] = compute_N2_from_profile(dens_500_pa[:,lat_i], depth)
#%% Plotting

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lat_NZ, depth, N2_lat_500_nz - N2_lat_1_nz, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax1.set_title('a) NZ convective region', fontsize=14)
ax1.set_xlim(-80,-30)
ax1.set_ylim(2000, 0)
cs = ax1.contour(lat_NZ, depth, dens_1_nz, levels = [1027.0], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lat_NZ, depth, dens_500_nz, levels = [1027.0], linestyles = '--', colors='black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lat_NZ, depth, dens_1_nz, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lat_NZ, depth, dens_500_nz, levels = [1027.5], linestyles = '--', colors='black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lat_NZ, depth, dens_1_nz, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lat_NZ, depth, dens_500_nz, levels = [1027.7], linestyles = '--', colors='black')
ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax1.vlines(x=-70, ymin=2000, ymax=0, color='grey')
ax1.vlines(x=-62, ymin=2000, ymax=0, color='grey')

CS2 = ax2.contourf(lat, depth, N2_lat_500_wgkp - N2_lat_1_wgkp, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax2.set_title('b) WGKP convective region', fontsize=14)
ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
cs = ax2.contour(lat, depth, dens_1_wgkp, levels = [1027.0], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_500_wgkp, levels = [1027.0], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1_wgkp, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_500_wgkp, levels = [1027.5], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_1_wgkp, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lat, depth, dens_500_wgkp, levels = [1027.7], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lat, depth, N2_lat_500_au - N2_lat_1_au, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax3.set_title('c) AU convective region', fontsize=14)
ax3.set_xlim(-80,-30)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lat, depth, dens_1_au, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_500_au, levels = [1027.0], linestyles = '--', colors='black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1_au, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_500_au, levels = [1027.5], linestyles = '--', colors='black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_1_au, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lat, depth, dens_500_au, levels = [1027.7], linestyles = '--', colors='black')
ax3.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax3.vlines(x=-70, ymin=2000, ymax=0, color='grey')
ax3.vlines(x=-40, ymin=2000, ymax=0, color='grey')

CS2 = ax4.contourf(lat, depth, N2_lat_500_pa - N2_lat_1_pa, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax4.set_title('d) PA convective region', fontsize=14)
ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
cs = ax4.contour(lat, depth, dens_1_pa, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500_pa, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1_pa, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500_pa, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_1_pa, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lat, depth, dens_500_pa, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Latitude [$^\circ$S]', fontsize=12)
ax4.vlines(x=-70, ymin=2000, ymax=0, color='grey')
ax4.vlines(x=-40, ymin=2000, ymax=0, color='grey')

plt.suptitle('N$^2$ anomalies ((500 - 600) minus (1 - 100))', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Buoyancy_frequency_convective_regions_zonal.pdf')

# %% Meridional transects of N2 anomalies

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_meridional_averaged_-70N--50N_AU_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_au         = fh.variables['lon'][:]    #longitude
temp_mer_1_au      = fh.variables['TEMP'][:]   #temperature
salt_mer_1_au      = fh.variables['SALT'][:]   #salinity
dens_mer_1_au      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_meridional_averaged_-70N--50N_AU_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_au         = fh.variables['lon'][:]    #longitude
temp_mer_500_au      = fh.variables['TEMP'][:]   #temperature
salt_mer_500_au      = fh.variables['SALT'][:]   #salinity
dens_mer_500_au      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_meridional_averaged_-70N--50N_PA_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_pa         = fh.variables['lon'][:]    #longitude
temp_mer_1_pa      = fh.variables['TEMP'][:]   #temperature
salt_mer_1_pa      = fh.variables['SALT'][:]   #salinity
dens_mer_1_pa      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_meridional_averaged_-70N--50N_PA_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_pa         = fh.variables['lon'][:]    #longitude
temp_mer_500_pa      = fh.variables['TEMP'][:]   #temperature
salt_mer_500_pa      = fh.variables['SALT'][:]   #salinity
dens_mer_500_pa      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_meridional_averaged_-70N--60N_NZ_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_nz         = fh.variables['lon'][:]    #longitude
temp_mer_1_nz      = fh.variables['TEMP'][:]   #temperature
salt_mer_1_nz      = fh.variables['SALT'][:]   #salinity
dens_mer_1_nz      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_meridional_averaged_-70N--60N_NZ_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_nz         = fh.variables['lon'][:]    #longitude
temp_mer_500_nz      = fh.variables['TEMP'][:]   #temperature
salt_mer_500_nz      = fh.variables['SALT'][:]   #salinity
dens_mer_500_nz      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_1-100_meridional_averaged_-80N--50N_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_wgkp         = fh.variables['lon'][:]    #longitude
temp_mer_1_wgkp      = fh.variables['TEMP'][:]   #temperature
salt_mer_1_wgkp      = fh.variables['SALT'][:]   #salinity
dens_mer_1_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()

fh = netcdf.Dataset(directory_data + 'TEMP_SALT_DENS_year_500-600_meridional_averaged_-80N--50N_WGKP_SO.nc','r')

depth       = fh.variables['depth'][:]  #depth
lon_wgkp         = fh.variables['lon'][:]    #longitude
temp_mer_500_wgkp      = fh.variables['TEMP'][:]   #temperature
salt_mer_500_wgkp      = fh.variables['SALT'][:]   #salinity
dens_mer_500_wgkp      = fh.variables['PD'][:]     #potential density

fh.close()

# %% Compute N2 for each profile

n0_lon_1_wgkp = ma.masked_all((len(depth), len(lon_wgkp)))
n0_lon_500_wgkp = ma.masked_all((len(depth), len(lon_wgkp)))
n0_lon_1_nz = ma.masked_all((len(depth), len(lon_nz)))
n0_lon_500_nz = ma.masked_all((len(depth), len(lon_nz)))
n0_lon_1_au = ma.masked_all((len(depth), len(lon_au)))
n0_lon_500_au = ma.masked_all((len(depth), len(lon_au)))
n0_lon_1_pa = ma.masked_all((len(depth), len(lon_pa)))
n0_lon_500_pa = ma.masked_all((len(depth), len(lon_pa)))

N2_lon_1_wgkp = ma.masked_all((len(depth), len(lon_wgkp)))
N2_lon_500_wgkp = ma.masked_all((len(depth), len(lon_wgkp)))
N2_lon_1_nz = ma.masked_all((len(depth), len(lon_nz)))
N2_lon_500_nz = ma.masked_all((len(depth), len(lon_nz)))
N2_lon_1_au = ma.masked_all((len(depth), len(lon_au)))
N2_lon_500_au = ma.masked_all((len(depth), len(lon_au)))
N2_lon_1_pa = ma.masked_all((len(depth), len(lon_pa)))
N2_lon_500_pa = ma.masked_all((len(depth), len(lon_pa)))

for lon_i in range(len(lon_nz)):
    n0_lon_1_nz[:, lon_i], N2_lon_1_nz[:,lon_i] = compute_N2_from_profile(dens_mer_1_nz[:,lon_i], depth)
    n0_lon_500_nz[:, lon_i], N2_lon_500_nz[:,lon_i] = compute_N2_from_profile(dens_mer_500_nz[:,lon_i], depth)

for lon_i in range(len(lon_au)):
    n0_lon_1_au[:, lon_i], N2_lon_1_au[:,lon_i] = compute_N2_from_profile(dens_mer_1_au[:,lon_i], depth)
    n0_lon_500_au[:, lon_i], N2_lon_500_au[:,lon_i] = compute_N2_from_profile(dens_mer_500_au[:,lon_i], depth)

for lon_i in range(len(lon_pa)):
    n0_lon_1_pa[:, lon_i], N2_lon_1_pa[:,lon_i] = compute_N2_from_profile(dens_mer_1_pa[:,lon_i], depth)
    n0_lon_500_pa[:, lon_i], N2_lon_500_pa[:,lon_i] = compute_N2_from_profile(dens_mer_500_pa[:,lon_i], depth)

for lon_i in range(len(lon_wgkp)):
    n0_lon_1_wgkp[:, lon_i], N2_lon_1_wgkp[:,lon_i] = compute_N2_from_profile(dens_mer_1_wgkp[:,lon_i], depth)
    n0_lon_500_wgkp[:, lon_i], N2_lon_500_wgkp[:,lon_i] = compute_N2_from_profile(dens_mer_500_wgkp[:,lon_i], depth)

# %%

fig, axs = plt.subplots(2, 2, figsize=(13, 8))

ax1, ax2, ax3, ax4 = axs.flatten()

CS2 = ax1.contourf(lon_nz, depth, N2_lon_500_nz - N2_lon_1_nz, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax1.set_title('a) NZ convective region', fontsize=14)
ax1.set_xlim(160,190)
ax1.set_ylim(2000, 0)
#cs = ax1.contour(lon_nz, depth, dens_mer_1_nz, levels = [1027.0], colors = 'black')
#ax1.clabel(cs, inline=True, manual=[(155, 1000)], fontsize=10, fmt='%.1f')
#cs = ax1.contour(lon_nz, depth, dens_mer_500_nz, levels = [1027.0], linestyles = '--', colors='black')
#ax1.clabel(cs, inline=True, manual=[(-60, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lon_nz, depth, dens_mer_1_nz, levels = [1027.5], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(155, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lon_nz, depth, dens_mer_500_nz, levels = [1027.5], linestyles = '--', colors='black')
ax1.clabel(cs, inline=True, manual=[(155, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lon_nz, depth, dens_mer_1_nz, levels = [1027.7], colors = 'black')
ax1.clabel(cs, inline=True, manual=[(170, 1000)], fontsize=10, fmt='%.1f')
cs = ax1.contour(lon_nz, depth, dens_mer_500_nz, levels = [1027.7], linestyles = '--', colors='black')
ax1.clabel(cs, inline=True, manual=[(170, 1000)], fontsize=10, fmt='%.1f')
ax1.set_ylabel('Depth [m]', fontsize=12)
ax1.set_xlabel('Longitude [$^\circ$E]', fontsize=12)

CS2 = ax2.contourf(lon_wgkp, depth, N2_lon_500_wgkp - N2_lon_1_wgkp, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax2.set_title('b) WGKP convective region', fontsize=14)
#ax2.set_xlim(-80,-30)
ax2.set_ylim(2000, 0)
#cs = ax2.contour(lon_wgkp, depth, dens_mer_1_wgkp, levels = [1027.0], colors = 'black')
#ax2.clabel(cs, inline=True, manual=[(130, 1000)], fontsize=10, fmt='%.1f')
#cs = ax2.contour(lon_wgkp, depth, dens_mer_500_wgkp, levels = [1027.0], linestyles = '--', colors='black')
#ax2.clabel(cs, inline=True, manual=[(130, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lon_wgkp, depth, dens_mer_1_wgkp, levels = [1027.5], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lon_wgkp, depth, dens_mer_500_wgkp, levels = [1027.5], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(60, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lon_wgkp, depth, dens_mer_1_wgkp, levels = [1027.7], colors = 'black')
ax2.clabel(cs, inline=True, manual=[(20, 1000)], fontsize=10, fmt='%.1f')
cs = ax2.contour(lon_wgkp, depth, dens_mer_500_wgkp, levels = [1027.7], linestyles = '--', colors='black')
ax2.clabel(cs, inline=True, manual=[(20, 1000)], fontsize=10, fmt='%.1f')
ax2.set_ylabel('Depth [m]', fontsize=12)
ax2.set_xlabel('Longitude [$^\circ$E]', fontsize=12)
#ax2.vlines(x=-80, ymin=2000, ymax=0, color='grey')
#ax2.vlines(x=-50, ymin=2000, ymax=0, color='grey')

CS2 = ax3.contourf(lon_au, depth, N2_lon_500_au - N2_lon_1_au, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax3.set_title('c) AU convective region', fontsize=14)
ax3.set_ylim(2000, 0)
cs = ax3.contour(lon_au, depth, dens_mer_1_au, levels = [1027.0], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(130, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lon_au, depth, dens_mer_500_au, levels = [1027.0], linestyles = '--', colors='black')
#ax3.clabel(cs, inline=True, manual=[(130, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lon_au, depth, dens_mer_1_au, levels = [1027.5], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(120, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lon_au, depth, dens_mer_500_au, levels = [1027.5], linestyles = '--', colors='black')
ax3.clabel(cs, inline=True, manual=[(120, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lon_au, depth, dens_mer_1_au, levels = [1027.7], colors = 'black')
ax3.clabel(cs, inline=True, manual=[(120, 1000)], fontsize=10, fmt='%.1f')
cs = ax3.contour(lon_au, depth, dens_mer_500_au, levels = [1027.7], linestyles = '--', colors='black')
ax3.clabel(cs, inline=True, manual=[(120, 1000)], fontsize=10, fmt='%.1f')
ax3.set_ylabel('Depth [m]', fontsize=12)
ax3.set_xlabel('Longitude [$^\circ$E]', fontsize=12)

CS2 = ax4.contourf(lon_pa, depth, N2_lon_500_pa - N2_lon_1_pa, levels= np.linspace(-4e-6, 4e-6, 21), extend = 'both', cmap = 'RdBu_r')

#plt.legend()
cbar	= colorbar(CS2, ticks = np.arange(-4e-6, 4e-6 + 1e-7, 1e-6))
cbar.set_label(r'$\Delta N²$ [s⁻²]', fontsize = 11)
ax4.set_title('d) PA convective region', fontsize=14)
#ax4.set_xlim(-80,-30)
ax4.set_ylim(2000, 0)
cs = ax4.contour(lon_pa, depth, dens_mer_1_pa, levels = [1027.0], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-80, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lon_pa, depth, dens_mer_500_pa, levels = [1027.0], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-80, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lon_pa, depth, dens_mer_1_pa, levels = [1027.5], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-80, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lon_pa, depth, dens_mer_500_pa, levels = [1027.5], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-80, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lon_pa, depth, dens_mer_1_pa, levels = [1027.7], colors = 'black')
ax4.clabel(cs, inline=True, manual=[(-80, 1000)], fontsize=10, fmt='%.1f')
cs = ax4.contour(lon_pa, depth, dens_mer_500_pa, levels = [1027.7], linestyles = '--', colors='black')
ax4.clabel(cs, inline=True, manual=[(-80, 1000)], fontsize=10, fmt='%.1f')
ax4.set_ylabel('Depth [m]', fontsize=12)
ax4.set_xlabel('Longitude [$^\circ$E]', fontsize=12)

plt.suptitle('Zonal N$^2$ anomalies ((500 - 600) minus (1 - 100))', fontsize=16)
plt.tight_layout()

fig.savefig(directory_figures +'Buoyancy_frequency_convective_regions_meridional.pdf')

 # %%
