#Program determines kinetic energy terms of the Lorenz Cycle of the SOM used for final computation

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


def ReadinData(filename, layer_avail = False, volume_norm = False):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#First get the u-grid
	lon 		= fh.variables['ULONG'][400]				#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['ULAT'][:, 780]				#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	dx		= fh.variables['DXU'][:]  / 100.0			#Zonal grid cell length (m)
	depth  	 	= fh.variables['z_t'][:]  / 100.0			#Depth (m)
	depth_grid 	= fh.variables['HU'][:] / 100.0				#Depth at u-grid (m)
	layer		= fh.variables['dz'][:] / 100.0				#Layer thickness (m)
	depth_top	= fh.variables['z_w_top'][:] / 100.0			#Top of grid cell
	area		= fh.variables['UAREA'][:] / 10000			#UAREA (m^2)

	fh.close()
	
	#print(lon)
	
	#plt.figure()
	#plt.plot(lon)
	#plt.show()
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0
	
	#plt.figure()
	#plt.plot(lon)
	#plt.show()

	#print(lon)
	#print(np.shape(lon))
	#print(np.min(lon))

	#Select region (SO30 region or WGKP region)
	#lon_min_index	= (fabs(lon - -35)).argmin()
	#lon_max_index	= (fabs(lon - 80)).argmin() + 1
	#lat_min_index	= (fabs(lat - -90)).argmin()
	#lat_max_index	= (fabs(lat - -50)).argmin() + 1
	
	#lon_min_index	= (fabs(lon - -35)).argmin()
	#lon_max_index	= (fabs(lon - -34)).argmin() + 1
	lat_min_index	= (fabs(lat - -90)).argmin()
	lat_max_index	= (fabs(lat - -30)).argmin() + 1
	
	#print(lat_min_index)
	#print(lat_max_index)	
		
	#lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, :]
	area		= area[lat_min_index:lat_max_index, :]
	
	#print(lon)
	#print(lat)
	
	##sys.exit()
		
	fh = netcdf.Dataset(filename, 'r')

	u_vel 		= fh.variables['UVEL'][:, lat_min_index:lat_max_index, :] / 100.0	#Zonal velocity (m/s)
	v_vel 		= fh.variables['VVEL'][:, lat_min_index:lat_max_index, :] / 100.0	#Meridional velocity (m/s)
	
	fh.close()
	
	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		u_vel[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], u_vel[depth_i])
		v_vel[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], v_vel[depth_i])

	#------------------------------------------------------------------------------
	if layer_avail	== False:	
		#Determine the depth per grid cell (parcel bottom cells)
		layer_field	= ma.masked_all((len(depth), len(lat), len(lon)))

		for depth_i in range(len(layer)):
			#print(depth_i)

			#Mask all elements which are land and fill the layer field with the depth layer for each layer
			u_vel_depth		= u_vel[depth_i]
			v_vel_depth		= v_vel[depth_i]
			u_vel_depth		= ma.masked_where((u_vel_depth == 0) & (v_vel_depth == 0), u_vel_depth) #mask cells where both uvel and vvel are zero (then it is some kind of bug in the HR pop)
			layer_field[depth_i]	= layer[depth_i]
			layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = u_vel_depth.mask)

			#Determine where the layer needs to be adjusted, partial depth cells
			depth_diff		= np.sum(layer_field, axis = 0) - depth_grid

			#If the depth difference is negative (i.e. bottom is not reached), set to zero
			depth_diff		= ma.masked_where(depth_diff < 0, depth_diff)
			depth_diff		= depth_diff.filled(fill_value = 0.0)

			#Subtract the difference of the current layer with the difference
			layer_field[depth_i]	= layer_field[depth_i] - depth_diff

		#Get the total vertical extent for each layer
		volume			= layer_field * area
		volume			= ma.masked_where(volume <= 0, volume) #volume of gridcells
		
	return lat, lon, depth, volume, u_vel, v_vel
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#First SOM cycle (model year 63-114), second (324-378) or last SOM cycle (500-600)
year_start	= 500
year_end	= 600

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
lat, lon, depth, volume, u_vel, v_vel		= ReadinData(files[0])
time_year					= ma.masked_all(year_end-year_start+1)
u_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
v_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
uu_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
vv_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))

#Volume integrated energetics
TKE_int  = ma.masked_all((len(time_year) - window + 1))
MKE_int  = ma.masked_all((len(time_year) - window + 1))
EKE_int  = ma.masked_all((len(time_year) - window + 1))

#------------------------------------------------------------------------------
month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.]*window)
month_days	= month_days / np.sum(month_days)

#Fill the array's with the same dimensions
month_days_all	= ma.masked_all((window*12, len(lat), len(lon)))

for month_i in range(len(month_days)):
	month_days_all[month_i]		= month_days[month_i]
	
#------------------------------------------------------------------------------
	
rho0 = 1026 #global average density of sea water [kg/m^3]

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start
		
	for month_i in range(12):
		
		#print(year_i+300+year_start-1)
		
		#Get the monthly files
		filename 	= files[year_i*12 + month_i]

		#print(filename)
		lat, lon, depth, volume, u_vel_year, v_vel_year = ReadinData(filename)
		
		# Compute the rolling index in uu_vel_all for the current month
		pdw_index = (year_i % window) * 12 + month_i
		#print(f"Year: {year_i}, Month: {month_i}, PDW_all index: {pdw_index}")
		
		#Save monthly u*u in window
		uu_vel_all[pdw_index,:,:,:] = u_vel_year**2
		vv_vel_all[pdw_index,:,:,:] = v_vel_year**2
		
		u_vel_all[pdw_index,:,:,:] = u_vel_year
		v_vel_all[pdw_index,:,:,:] = v_vel_year
		
	#Start calculating energetics when year_i has at least 5 years
	if year_i >= window - 1:

		#5-year window available
		TKE_year	= ma.masked_all((len(depth), len(lat), len(lon)))	#Total kinetic energy
		MKE_year	= ma.masked_all((len(depth), len(lat), len(lon)))	#Mean kinetic energy	
		EKE_year	= ma.masked_all((len(depth), len(lat), len(lon)))	#Eddy (or residual) kinetic energy
		
		for depth_i in range(len(depth)):
			TKE_year[depth_i]	= 0.5 * rho0 * np.sum((uu_vel_all[:, depth_i] + vv_vel_all[:, depth_i])*month_days_all, axis = 0) 		
			MKE_year[depth_i]	= 0.5 * rho0 * (np.sum(u_vel_all[:, depth_i]*month_days_all, axis = 0)**2.0 + np.sum(v_vel_all[:, depth_i]*month_days_all, axis = 0)**2.0) 				
			EKE_year[depth_i]	= TKE_year[depth_i] - MKE_year[depth_i]										
		
		#Volume integrated
		TKE_int[year_i - window + 1] = np.sum(TKE_year*volume)
		MKE_int[year_i - window + 1] = np.sum(MKE_year*volume)
		EKE_int[year_i - window + 1] = np.sum(EKE_year*volume)

	#-----------------------------------------------------------------------------------------
	#-----------------------------------------------------------------------------------------
	#continue
	print('Data is written to file')
	fh = netcdf.Dataset(directory+'Ocean/KE_year_'+str(year_start)+'-'+str(year_end)+'_window_'+str(window)+'_volume_integrated_SO30_finalcode.nc', 'w')

	fh.createDimension('time', len(time_year) - window + 1)
	fh.createDimension('lat', len(lat))
	fh.createDimension('lon', len(lon))

	fh.createVariable('time', float, ('time'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('TKE', float, ('time'), zlib=True)
	fh.createVariable('MKE', float, ('time'), zlib=True)
	fh.createVariable('EKE', float, ('time'), zlib=True)

	fh.variables['time'].long_name 	= 'Starting year of window'
	fh.variables['lat'].long_name	= 'Latitudes for surface integration'
	fh.variables['lon'].long_name	= 'Longitudes for surface integration'
	fh.variables['TKE'].long_name 	= 'Volume integrated total kinetic energy'
	fh.variables['MKE'].long_name 	= 'Volume integrated mean kinetic energy'
	fh.variables['EKE'].long_name 	= 'Volume integrated eddy kinetic energy'

	fh.variables['time'].units 	= 'year'
	fh.variables['lat'].units 	= 'Degrees N'
	fh.variables['lon'].units 	= 'Degrees E'
	fh.variables['TKE'].units 	= 'J'
	fh.variables['MKE'].units 	= 'J'
	fh.variables['EKE'].units 	= 'J'

	#Writing data to correct variable
	fh.variables['time'][:] 	= time_year[:len(time_year) - window + 1]
	fh.variables['lat'][:] 		= lat
	fh.variables['lon'][:] 		= lon
	fh.variables['TKE'][:] 		= TKE_int
	fh.variables['MKE'][:] 		= MKE_int
	fh.variables['EKE'][:] 		= EKE_int

	fh.close()	
	
	
	
	
	
	
	
