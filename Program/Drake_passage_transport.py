#Program reads-in three dimensional velocity field
#The mass transport is determined at 26N (AMOC)

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

def ReadinData(filename, lon_index, lat_min_index, lat_max_index, depth_min_index, depth_max_index):

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	lon 		= fh.variables['ULONG'][lat_min_index:lat_max_index, lon_index]			#Longitude
	lat 		= fh.variables['ULAT'][lat_min_index:lat_max_index, lon_index]			#Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] / 100.0			#Depth (m)
	depth_u 	= fh.variables['HU'][lat_min_index:lat_max_index, lon_index] / 100.0		#Depth at u-grid (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] / 100.0			#Layer thickness (m)
	grid_y		= fh.variables['DYU'][lat_min_index:lat_max_index, lon_index] / 100.0		#Meridional grid cell length (m)
	grid_yt		= fh.variables['DYT'][lat_min_index:lat_max_index, lon_index] / 100.0
	depth_top	= fh.variables['z_w_top'][:] / 100.0						#Top of grid cell
		
	fh.close()
	
	plt.figure()
	plt.plot(lat, grid_y)
	plt.plot(lat, grid_yt)
	plt.show()
	
	sys.exit()
	
	fh = netcdf.Dataset(filename, 'r')

	u_vel 		= fh.variables['UVEL'][depth_min_index:depth_max_index, lat_min_index:lat_max_index, lon_index] / 100.0	#Zonal velocity (m/s)

	fh.close()

	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		u_vel[depth_i]	= ma.masked_where(depth_u <= depth_top[depth_i], u_vel[depth_i])

	return lon, lat, depth, layer, depth_u, grid_y, u_vel
			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

num_year	= 600

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset('/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc', 'r')

lon 		= fh.variables['ULONG'][:]		#Longitude
lat 		= fh.variables['ULAT'][:]		#Latitude  
depth   	= fh.variables['z_t'][:] / 100.0	#Depth (m)
depth_grid	= fh.variables['HU'][:] / 100.0		#Depth at u-grid (m)
	
fh.close()

depth_min 	= 0
depth_max	= 6000

#Get the dimensions of depth, latitude and longitude
#The lon-lat indices are retained from the CESM Diagnostics page
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lon_index		= 436
lat_min_index 		= 287
lat_max_index 		= 522
depth_grid		= depth_grid[lat_min_index:lat_max_index, lon_index]

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, depth, layer, depth_u, grid_y, u_vel 	= ReadinData('/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/tavg/t.t0.1_42l_nccs01.030002.nc', lon_index, lat_min_index, lat_max_index, depth_min_index, depth_max_index)
layer_field				= ma.masked_all((len(depth), len(lat)))

print(lon)
print(lat)

sys.exit()

for depth_i in range(len(depth)):
	#Mask all elements which are land and fill the layer field with the depth layer for each layer
	layer_field[depth_i]	= layer[depth_i]
	layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = u_vel[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff		= np.sum(layer_field, axis = 0) - depth_u

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff		= ma.masked_where(depth_diff < 0, depth_diff)
	depth_diff		= depth_diff.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field[depth_i]	= layer_field[depth_i] - depth_diff
		
#-----------------------------------------------------------------------------------------

year_start	= 0
year_end	= 600

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

#Define empty array's
time_year	= ma.masked_all(num_year)
transport_all	= ma.masked_all(len(time_year))

for year_i in range(num_year):
	print(year_i)
	time_year[year_i] 	= year_i + 1
	u_vel	 		= ma.masked_all((12, len(depth), len(lat)))
	
	for month_i in range(12):

		#Get the monthly files
		filename 	= files[year_i*12 + month_i]

		print(filename)
		
		lon, lat, depth, layer, depth_u, grid_y, u_vel_month 	= ReadinData(filename, lon_index, lat_min_index, lat_max_index, depth_min_index, depth_max_index)

		u_vel[month_i]	= u_vel_month

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	u_vel	= np.sum(u_vel * month_days_all, axis = 0)

	#------------------------------------------------------------------------------

	#Determine the meridional transport
	transport	= u_vel * layer_field * grid_y

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[year_i]	= np.sum(transport) / 1000000.0
	
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'/Ocean/Drake_Passage_transport.nc', 'w')

fh.createDimension('time', len(time_year))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)

fh.variables['Transport'].longname 	= 'Volume transport'

fh.variables['time'].units 		= 'Year'
fh.variables['Transport'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time_year
fh.variables['Transport'][:]    	= transport_all

fh.close()
