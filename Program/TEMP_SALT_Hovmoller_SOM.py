#Program determines the MOV index

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

def ReadinData(filename, depth_min, depth_max, layer_avail = False, volume_norm = False):
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

	#Select region
	depth_min_index	= (fabs(depth_min - depth)).argmin()
	depth_max_index	= (fabs(depth_max - depth)).argmin() + 1
	depth		= depth[depth_min_index:depth_max_index]
	layer		= layer[depth_min_index:depth_max_index]
	depth_top	= depth_top[depth_min_index:depth_max_index]
	lon_min_index	= (fabs(lon - -55)).argmin()
	lon_max_index	= (fabs(lon - -5)).argmin() + 1
	lat_min_index	= (fabs(lat - -75)).argmin()
	lat_max_index	= (fabs(lat - 10)).argmin() + 1	
	
	lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
		
	fh = netcdf.Dataset(filename, 'r')

	temp 		= fh.variables['TEMP'][depth_min_index:depth_max_index, lat_min_index:lat_max_index, lon_min_index:lon_max_index]	#Potential temperature (deg C)
	salt 		= fh.variables['SALT'][depth_min_index:depth_max_index, lat_min_index:lat_max_index, lon_min_index:lon_max_index] * 1000. #Salinity (g /kg)
	
	fh.close()

	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		temp[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], temp[depth_i])
		salt[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], salt[depth_i])
			
	#------------------------------------------------------------------------------
	if layer_avail	== False:
	
		fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

		depth   	= fh.variables['z_t'][:depth_max_index]  / 100.0	#Depth (m)
		layer		= fh.variables['dz'][:depth_max_index] / 100.0		#Layer thickness (m)
		depth_top	= fh.variables['z_w_top'][:depth_max_index] / 100.0	#Top of grid cell
		
		fh.close()
		
		fh = netcdf.Dataset(filename, 'r')

		temp 		= fh.variables['TEMP'][:depth_max_index, lat_min_index:lat_max_index, lon_min_index:lon_max_index]	#Potential temperature (deg C)

		fh.close()

		for depth_i in range(len(depth)):
			#Mask all the field at the topography
			temp[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], temp[depth_i])
	
		#Determine the depth per grid cell (parcel bottom cells)
		layer_field	= ma.masked_all((len(depth), len(lat), len(lon)))

		for depth_i in range(len(layer)):
			print(depth_i)

			#Mask all elements which are land and fill the layer field with the depth layer for each layer
			temp_depth		= temp[depth_i]
			layer_field[depth_i]	= layer[depth_i]
			layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = temp_depth.mask)

			#Determine where the layer needs to be adjusted, partial depth cells
			depth_diff		= np.sum(layer_field, axis = 0) - depth_grid

			#If the depth difference is negative (i.e. bottom is not reached), set to zero
			depth_diff		= ma.masked_where(depth_diff < 0, depth_diff)
			depth_diff		= depth_diff.filled(fill_value = 0.0)

			#Subtract the difference of the current layer with the difference
			layer_field[depth_i]	= layer_field[depth_i] - depth_diff
			
		

		#------------------------------------------------------------------------------
		#Empty array for the depth indices and for the weighted layers
		layer_field_regions	= ma.masked_all((len(depth), len(lat), len(lon)))

		#Get the depth indices for each region
		depth_min_index = np.where(depth_top <= depth_min)[0][-1]
		depth_max_index	= np.where(depth_top <= depth_max)[0][-1] + 1

		depth_sections	= np.zeros(depth_max_index - depth_min_index)

		for depth_i in range(depth_min_index, depth_max_index):
			#Loop over the depth interval

			if depth_i == depth_min_index:
				#First layer, get the depth difference with respect to top
				depth_diff		= depth_min - depth_top[depth_i]

				depth_sections[depth_i - depth_min_index] = layer[depth_i] - depth_diff

			elif depth_i == depth_max_index - 1:
				#Last layer, get the depth difference with respect to top
				depth_diff	= depth_max -  depth_top[depth_i]

				depth_sections[depth_i - depth_min_index] = depth_diff

			else:
				#For all the other layers
				depth_sections[depth_i - depth_min_index] = layer[depth_i]

		#------------------------------------------------------------------------------
		for depth_i in range(depth_min_index, depth_max_index):
			#Now only get the relevant region and depth interval
			layer_region	= layer_field[depth_i] 

			if depth_i == depth_min_index:
				#Upper layer
				#Set the regions which are above the depth section to masked elements
				layer_region	= ma.masked_where(layer[depth_i] - layer_region > depth_sections[depth_i - depth_min_index], layer_region)

				#Set the remaining depth to bottom as layer thickness
				layer_region	= depth_sections[depth_i - depth_min_index] - (layer[depth_i] - layer_region)

			else:
				#The other layers
				#Get the layers where the bottom of the layer is not reached
				index_below	= np.where(layer_region < depth_sections[depth_i - depth_min_index])
				index_all	= np.where(layer_region >= depth_sections[depth_i - depth_min_index])

				layer_region[index_below]	= layer_region[index_below]
				layer_region[index_all]		= layer_region[index_all] * 0.0 + depth_sections[depth_i - depth_min_index]

			#Save the depth layer per grid cell
			layer_field_regions[depth_i]	= layer_region
		
		#Remove the depths which are not needed
		depth_min_index		= (fabs(depth_min - depth)).argmin()		
		layer_field_regions	= layer_field_regions[depth_min_index:]
		temp			= temp[depth_min_index:]
		depth			= depth[depth_min_index:]
		layer			= layer[depth_min_index:]
		depth_top		= depth_top[depth_min_index:]

		#Get the total vertical extent for each layer
		volume			= layer_field_regions * area
		volume_norm		= ma.masked_all(np.shape(volume))
		
		for lat_i in range(len(lat)):	
			#Normalise the field for each latitude
			volume_norm[:, lat_i]	= volume[:, lat_i] / np.sum(volume[:, lat_i])
			
	#Take the volume-averaged zonal mean for temperature
	temp		= np.sum(temp * volume_norm, axis = (0, 2))
	salt		= np.sum(salt * volume_norm, axis = (0, 2))
	
	return lat, volume_norm, temp, salt
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 300
depth_max	= 700

year_start	= 1
year_end	= 600

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

files	= files[:600*12]

#-----------------------------------------------------------------------------------------
lat, volume_norm, temp, salt	= ReadinData(files[0], depth_min, depth_max)
time_year			= ma.masked_all(year_end-year_start+1)
temp_all			= ma.masked_all((len(time_year), len(lat)))
salt_all			= ma.masked_all((len(time_year), len(lat)))

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	temp_year 		= ma.masked_all((12, len(lat)))
	salt_year 		= ma.masked_all((12, len(lat)))
		
	for month_i in range(12):

		#Get the monthly files 
    filename 	= files[year_i*12 + month_i]

		print(filename)
		lat, volume_norm, temp_year[month_i], salt_year[month_i] = ReadinData(filename, depth_min, depth_max, True, volume_norm)

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(lat)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	temp_all[year_i]	= np.sum(temp_year * month_days_all, axis = 0)
	salt_all[year_i]	= np.sum(salt_year * month_days_all, axis = 0)
	
#temp_plot	= temp_all - np.mean(temp_all, axis = 0)
	
#contourf(time_year, lat, salt_all.transpose())
#show()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/TEMP_SALT_Hovmoller_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'w')

fh.createDimension('time', len(time_year))
fh.createDimension('lat', len(lat))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('TEMP', float, ('time', 'lat'), zlib=True)
fh.createVariable('SALT', float, ('time', 'lat'), zlib=True)

fh.variables['lat'].long_name 		= 'Array of t-latitudes'
fh.variables['TEMP'].long_name 		= 'Potential temperature'
fh.variables['SALT'].long_name 		= 'Salinity'

fh.variables['time'].units 		= 'year'
fh.variables['lat'].units 		= 'Degrees N'
fh.variables['TEMP'].units 		= 'deg C'
fh.variables['SALT'].units 		= 'g / kg'

#Writing data to correct variable
fh.variables['time'][:] 		= time_year
fh.variables['lat'][:] 			= lat
fh.variables['TEMP'][:] 		= temp_all
fh.variables['SALT'][:] 		= salt_all

fh.close()
