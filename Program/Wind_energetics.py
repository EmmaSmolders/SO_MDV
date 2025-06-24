#Program determines wind energy input for SO30 region

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory_data	= '/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/tavg/'
directory	= '/home/smolders/HR_POP/Data/'

def ReadinData(filename):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#Grid is rectangular up to 28N (TAUX and TAUY are on U-grid)
	lon 		= fh.variables['ULONG'][400]		#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['ULAT'][:, 780]		#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	depth_u		= fh.variables['HU'][:] / 100		#Depth of bathymetry
	area		= fh.variables['UAREA'][:] / 10000	#TAREA (m^2)

	fh.close()
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0
	
	#Select region (WGKP region)
	#lon_min_index	= (fabs(lon - -35)).argmin()
	#lon_max_index	= (fabs(lon - 80)).argmin() + 1
	#lat_min_index	= (fabs(lat - -90)).argmin()
	#lat_max_index	= (fabs(lat - -50)).argmin() + 1
	
	#lon		= lon[lon_min_index:lon_max_index]
	#lat		= lat[lat_min_index:lat_max_index]
	#depth_u		= depth_u[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	#area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]

	#Select region (SO30)
	lat_min_index	= (fabs(lat - -90)).argmin()
	lat_max_index	= (fabs(lat - -30)).argmin() + 1	
	
	lat		= lat[lat_min_index:lat_max_index]
	depth_u	= depth_u[lat_min_index:lat_max_index, :]
	area		= area[lat_min_index:lat_max_index, :]
	
	fh = netcdf.Dataset(filename, 'r')

	#Get the relevant region
	tau_x 		= fh.variables['TAUX'][lat_min_index:lat_max_index, :]		#Horizontal windstress (dyne/cm^2 = g/(cm * s^2))
	tau_y 		= fh.variables['TAUY'][lat_min_index:lat_max_index, :]		#Meridional windstress (dyne/cm^2)
	u_vel 		= fh.variables['UVEL'][0, lat_min_index:lat_max_index, :] 		#Surface zonal velocity (cm/s)
	v_vel 		= fh.variables['VVEL'][0, lat_min_index:lat_max_index, :] 		#Surface meridional velocity (cm/s)

	fh.close()
	
	print(np.shape(depth_u))
	print(np.shape(tau_x))
	
	#Remove land mask
	tau_x		= ma.masked_where(depth_u <= 0, tau_x)
	tau_y		= ma.masked_where(depth_u <= 0, tau_y)
	u_vel		= ma.masked_where(depth_u <= 0, u_vel)
	v_vel		= ma.masked_where(depth_u <= 0, v_vel)	

	return lat, lon, area, tau_x, tau_y, u_vel, v_vel
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#Second SOM cycle (model year 63-114), second (324-378) or last SOM cycle (500-600)
year_start	= 324
year_end	= 378

#Moving average window
window = 5

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

#note that january 300 does not exist and that february 300 does not have PD data. 
files	= files[(year_start-1)*12:year_end*12]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------

lat, lon, area, tau_x, tau_y, u_vel, v_vel	= ReadinData(files[0])
time_year					= np.arange(year_start, year_end+1)
tau_x_all					= ma.masked_all((len(time_year), len(lat), len(lon)))
tau_y_all					= ma.masked_all((len(time_year), len(lat), len(lon)))
u_vel_all					= ma.masked_all((len(time_year), len(lat), len(lon)))
v_vel_all					= ma.masked_all((len(time_year), len(lat), len(lon)))
u_vel_tau_x_all					= ma.masked_all((len(time_year)*12, len(lat), len(lon)))
v_vel_tau_y_all					= ma.masked_all((len(time_year)*12, len(lat), len(lon)))

#-----------------------------------------------------------------------------------------

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	tau_x_year 		= ma.masked_all((12, len(lat), len(lon)))
	tau_y_year 		= ma.masked_all((12, len(lat), len(lon)))
	u_vel_year 		= ma.masked_all((12, len(lat), len(lon)))
	v_vel_year 		= ma.masked_all((12, len(lat), len(lon)))
	
	for month_i in range(12):
		
		print(year_i+300+year_start-1)
		
		#Get the monthly files (data from last 90 years for forward hosing is in subfolder)
		filename 	= files[year_i*12 + month_i]
		print(filename)

		print(filename)
		lat, lon, area, tau_x_year[month_i], tau_y_year[month_i], u_vel_year[month_i], v_vel_year[month_i] = ReadinData(filename)
		
		u_vel_tau_x_all[year_i*12 + month_i] = u_vel_year[month_i]*tau_x_year[month_i]
		v_vel_tau_y_all[year_i*12 + month_i] = v_vel_year[month_i]*tau_y_year[month_i]

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(lat), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	tau_x_all[year_i]	= np.sum(tau_x_year * month_days_all, axis = 0)
	tau_y_all[year_i]	= np.sum(tau_y_year * month_days_all, axis = 0)
	u_vel_all[year_i]	= np.sum(u_vel_year * month_days_all, axis = 0)
	v_vel_all[year_i]	= np.sum(v_vel_year * month_days_all, axis = 0)

#-----------------------------------------------------------------------------------------

gKm  		= ma.masked_all((len(time_year) - window + 1, len(lat), len(lon)))
gKe  		= ma.masked_all((len(time_year) - window + 1, len(lat), len(lon)))
gKtotal  	= ma.masked_all((len(time_year) - window + 1, len(lat), len(lon)))

for time_i in range(len(time_year) - window + 1):
	#Determine mean and eddy kinetic energy generation
	gKm[time_i, :, :]	= (np.mean(u_vel_all[time_i:time_i+window, :, :], axis = 0)*np.mean(tau_x_all[time_i:time_i+window, :, :], axis = 0) + np.mean(v_vel_all[time_i:time_i+window, :, :], axis = 0)*np.mean(tau_y_all[time_i:time_i+window, :, :], axis = 0))*1.0e-03	#Mean kinetic energy generation [kg/s^3]
	gKtotal[time_i, :, :]	= np.mean(u_vel_tau_x_all[time_i*12:(time_i+window)*12, :, :] + v_vel_tau_y_all[time_i*12:(time_i+window)*12, :, :], axis=0)*1.0e-03
	gKe[time_i, :, :]	= gKtotal[time_i, :, :] - gKm[time_i, :, :]	#Eddy kinetic energy generation [kg/s^3]

#Surface integrated energetics
gKm_int  	= ma.masked_all((len(time_year) - window + 1))
gKe_int  	= ma.masked_all((len(time_year) - window + 1))
gKtotal_int 	= ma.masked_all((len(time_year) - window + 1))

for time_i in range(len(time_year) - window + 1):
	gKm_int[time_i] 	= np.sum(gKm[time_i]*area) #[kg*m^2/s^3 = W]
	gKe_int[time_i] 	= np.sum(gKe[time_i]*area)
	gKtotal_int[time_i] 	= np.sum(gKtotal[time_i]*area)
	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/Generation_KE_year_'+str(year_start)+'-'+str(year_end)+'_window_'+str(window)+'_area_integrated_SO30.nc', 'w')

fh.createDimension('time', len(time_year) - window + 1)
fh.createDimension('lat', len(lat))
fh.createDimension('lon', len(lon))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('lon', float, ('lon'), zlib=True)
fh.createVariable('gKm', float, ('time'), zlib=True)
fh.createVariable('gKe', float, ('time'), zlib=True)
fh.createVariable('gKtotal', float, ('time'), zlib=True)

fh.variables['time'].long_name 		= 'Starting year of window'
fh.variables['lat'].long_name		= 'Latitudes for surface integration'
fh.variables['lon'].long_name		= 'Longitudes for surface integration'
fh.variables['gKm'].long_name 		= 'Area integrated mean kinetic energy generation'
fh.variables['gKe'].long_name 		= 'Area integrated eddy kinetic energy generation'
fh.variables['gKtotal'].long_name 	= 'Area integrated total kinetic energy generation'

fh.variables['time'].units 	= 'year'
fh.variables['gKe'].units 	= 'W'
fh.variables['gKm'].units 	= 'W'
fh.variables['gKtotal'].units 	= 'W'
fh.variables['lat'].units 	= 'Degrees N'
fh.variables['lon'].units 	= 'Degrees E'

#Writing data to correct variable
fh.variables['time'][:] 	= time_year[:len(time_year) - window + 1]
fh.variables['gKm'][:] 		= gKm_int
fh.variables['gKe'][:] 		= gKe_int
fh.variables['gKtotal'][:] 	= gKtotal_int
fh.variables['lat'][:] 		= lat
fh.variables['lon'][:] 		= lon

fh.close()	
	
	
	
	
