#Program determines area averaged temperature or potential density in WKGP region

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

	#Grid is rectangular up to 28N
	lon 		= fh.variables['TLONG'][400]		#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['TLAT'][:, 780]		#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	depth_grid	= fh.variables['HT'][:] / 100		#Depth of bathymetry
	area		= fh.variables['TAREA'][:] / 10000	#TAREA (m^2)
	depth   	= fh.variables['z_t'][:]  / 100.0	#Depth (m)
	layer		= fh.variables['dz'][:] / 100.0		#Layer thickness (m)
	depth_top	= fh.variables['z_w_top'][:] / 100.0	#Top of grid cell
	
	fh.close()
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0

	#Select region (WGKP region)
	lon_min_index	= (fabs(lon - -35)).argmin()
	lon_max_index	= (fabs(lon - 80)).argmin() + 1
	lat_min_index	= (fabs(lat - -90)).argmin()
	lat_max_index	= (fabs(lat - -50)).argmin() + 1
	
	lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
		
	fh = netcdf.Dataset(filename, 'r')

	#temp 		= fh.variables['TEMP'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index]		#Potential temperature (deg C)
	
	if 'PD' not in fh.variables:
		PD = fh.variables['RHO'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index] * 1000.  # In-situ density (kg/m^3)
	else:
		PD = fh.variables['PD'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index] * 1000.  # Potential density (kg/m^3)
	
	fh.close()
	
	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		#temp[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], temp[depth_i])
		PD[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], PD[depth_i])

	#------------------------------------------------------------------------------
	if layer_avail	== False:	
		#Determine the depth per grid cell (parcel bottom cells)
		layer_field	= ma.masked_all((len(depth), len(lat), len(lon)))

		for depth_i in range(len(layer)):
			#print(depth_i)

			#Mask all elements which are land and fill the layer field with the depth layer for each layer
			PD_depth		= PD[depth_i]
			PD_depth		= ma.masked_where(PD_depth <= 0, PD_depth)
			layer_field[depth_i]	= layer[depth_i]
			layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = PD_depth.mask)

			#Determine where the layer needs to be adjusted, partial depth cells
			depth_diff		= np.sum(layer_field, axis = 0) - depth_grid

			#If the depth difference is negative (i.e. bottom is not reached), set to zero
			depth_diff		= ma.masked_where(depth_diff < 0, depth_diff)
			depth_diff		= depth_diff.filled(fill_value = 0.0)

			#Subtract the difference of the current layer with the difference
			layer_field[depth_i]	= layer_field[depth_i] - depth_diff

		#Get the total vertical extent for each layer
		volume			= layer_field * area
		volume			= ma.masked_where(volume <= 0, volume)
		volume_norm		= ma.masked_all(np.shape(volume))
		
		for depth_i in range(len(depth)):
			for lat_i in range(len(lat)):	
				#Normalise the field for each latitude and depth layer
				volume_norm[depth_i, lat_i]	= volume[depth_i, lat_i] / np.sum(volume[depth_i, lat_i])
				
	#Take the volume-averaged zonal mean for temperature
	PD		= np.sum(PD * volume_norm, axis = 2)
	
	#Take also mean in latitudinal direction
	PD = np.nanmean(PD, axis=1)
	
	#print(np.min(dens))
	
	#CS	= contourf(lat, depth, dens, levels = np.arange(1024, 1030.1, 0.1), extend = 'both', cmap = 'Spectral_r')
	#colorbar(CS)
	#show()
	#sys.exit()
		
	return depth, volume_norm, PD
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#First year has no PD or RHO data
year_start	= 2
year_end	= 600

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

#Jan 300 is stored in 300-02
files	= files[(year_start-1)*12:year_end*12]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------
depth, volume_norm, PD	= ReadinData(files[0])
time_year		= ma.masked_all(year_end-year_start+1)
PD_all			= ma.masked_all((len(time_year), len(depth)))

print(np.shape(PD))

for year_i in range(len(time_year)):
	
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	PD_year 		= ma.masked_all((12, len(depth)))
		
	for month_i in range(12):

		#Get the monthly files
		filename 	= files[year_i*12 + month_i]

		print(filename)
		depth, volume_norm, PD_year[month_i] = ReadinData(filename)

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	PD_all[year_i]	= np.sum(PD_year * month_days_all, axis = 0)
	
	print(PD_all[year_i])

	print(PD_all)

	#-----------------------------------------------------------------------------------------
	#-----------------------------------------------------------------------------------------

	print('Data is written to file')
	fh = netcdf.Dataset(directory+'Ocean/PD_year_'+str(year_start)+'-'+str(year_end)+'_area_averaged_WGKP.nc', 'w')

	fh.createDimension('time', len(time_year))
	fh.createDimension('depth', len(depth))

	fh.createVariable('depth', float, ('depth'), zlib=True)
	fh.createVariable('time', float, ('time'), zlib=True)
	fh.createVariable('PD', float, ('time', 'depth'), zlib=True)

	fh.variables['time'].long_name 		= 'Time'
	fh.variables['depth'].long_name 	= 'Depth'
	fh.variables['PD'].long_name 		= 'Yearly, area averaged potential density'

	fh.variables['depth'].units 		= 'm'
	fh.variables['time'].units 		= 'year'
	fh.variables['PD'].units 		= 'kg/m3'

	#Writing data to correct variable
	fh.variables['time'][:] 		= time_year
	fh.variables['depth'][:] 		= depth
	fh.variables['PD'][:] 			= PD_all

	fh.close()	
	
	
