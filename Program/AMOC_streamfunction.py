#Program determines the AMOC strength

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

def ReadinData(filename, lat_index):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#First get the u-grid
	lon 		= fh.variables['ULONG'][lat_index, 289:1000]		#Longitude
	lat 		= fh.variables['ULAT'][lat_index, 289:1000]		#Latitude 
	dx		= fh.variables['DXU'][lat_index, 289:1000]  / 100.0	#Zonal grid cell length (m)
	depth  	 	= fh.variables['z_t'][:]  / 100.0			#Depth (m)
	depth_u 	= fh.variables['HU'][lat_index, 289:1000] / 100.0	#Depth at u-grid (m)
	layer		= fh.variables['dz'][:] / 100.0				#Layer thickness (m)
	depth_top	= fh.variables['z_w_top'][:] / 100.0			#Top of grid cell

	fh.close()

	#Get the dimensions of depth
	depth_min_index 	= (np.abs(depth_min - depth)).argmin()
	depth_max_index 	= (np.abs(depth_max - depth)).argmin() + 1
	depth			= depth[depth_min_index:depth_max_index]
	layer			= layer[depth_min_index:depth_max_index]
	depth_top		= depth_top[depth_min_index:depth_max_index]

	fh = netcdf.Dataset(filename, 'r')

	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_index, 289:1000] / 100.0	#Meridional velocity (m/s)

	fh.close()

	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		v_vel[depth_i]	= ma.masked_where(depth_u <= depth_top[depth_i], v_vel[depth_i])

	return lon, lat, depth, layer, depth_top, depth_u, dx, v_vel
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 1000

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

files	= files[:600*12]


#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-9:-3]	
	year  = int(date[0:4])
	month = int(date[4:6])

	time[year_i] = year + (month-1) / 12.0
#-----------------------------------------------------------------------------------------

#Lat index along 26N
lat_index		= 1450

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, depth, layer, depth_top, depth_u, grid_x, v_vel 	= ReadinData(files[0], lat_index)
layer_field							= ma.masked_all((len(depth), len(lon)))

print(np.shape(v_vel))
print(np.shape(lat))
print(np.shape(lon))

plt.plot(lon, lat)
show()

sys.exit()

for depth_i in range(len(layer)):

	#Mask all elements which are land and fill the layer field with the depth layer for each layer
	layer_field[depth_i]	= layer[depth_i]
	layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = v_vel[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff		= np.sum(layer_field, axis = 0) - depth_u

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff		= ma.masked_where(depth_diff < 0, depth_diff)
	depth_diff		= depth_diff.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field[depth_i]	= layer_field[depth_i] - depth_diff

	#------------------------------------------------------------------------------
	#Empty array for the depth indices and for the weighted layers
	layer_field_regions	= ma.masked_all((len(depth), len(lon)))

	#Get the depth indices for each region
	depth_min_index_2 = np.where(depth_top <= depth_min)[0][-1]
	depth_max_index_2 = np.where(depth_top <= depth_max)[0][-1] + 1

	depth_sections	= np.zeros(depth_max_index_2 - depth_min_index_2)

	for depth_i in range(depth_min_index_2, depth_max_index_2):
		#Loop over the depth interval

		if depth_i == depth_min_index_2:
			#First layer, get the depth difference with respect to top
			depth_diff		= depth_min - depth_top[depth_i]

			depth_sections[depth_i - depth_min_index_2] = layer[depth_i] - depth_diff

		elif depth_i == depth_max_index_2 - 1:
			#Last layer, get the depth difference with respect to top
			depth_diff	= depth_max -  depth_top[depth_i]

			depth_sections[depth_i - depth_min_index_2] = depth_diff

		else:
			#For all the other layers
			depth_sections[depth_i - depth_min_index_2] = layer[depth_i]

	#------------------------------------------------------------------------------
	for depth_i in range(depth_min_index_2, depth_max_index_2):
		#Now only get the relevant region and depth interval
		layer_region	= layer_field[depth_i] 

		if depth_i == depth_min_index_2:
			#Upper layer
			#Set the regions which are above the depth section to masked elements
			layer_region	= ma.masked_where(layer[depth_i] - layer_region > depth_sections[depth_i - depth_min_index_2], layer_region)

			#Set the remaining depth to bottom as layer thickness
			layer_region	= depth_sections[depth_i - depth_min_index_2] - (layer[depth_i] - layer_region)

		else:
			#The other layers
			#Get the layers where the bottom of the layer is not reached
			index_below	= np.where(layer_region < depth_sections[depth_i - depth_min_index_2])
			index_all	= np.where(layer_region >= depth_sections[depth_i - depth_min_index_2])

			layer_region[index_below]	= layer_region[index_below]
			layer_region[index_all]		= layer_region[index_all] * 0.0 + depth_sections[depth_i - depth_min_index_2]

		#Save the depth layer per grid cell
		layer_field_regions[depth_i]	= layer_region

#-----------------------------------------------------------------------------------------
time_year		= ma.masked_all(int(len(time)/12))
transport_all		= ma.masked_all(len(time_year))

for year_i in range(int(np.min(time)), int(np.min(time))+len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i - int(np.min(time))] = year_i - 300 + 1

	v_vel	=    ma.masked_all((12, len(depth), len(lon)))

	for month_i in range(12):

		#Get the monthly files
		filename	= directory_data+'t.t0.1_42l_nccs01.'+str(year_i).zfill(4)+str(month_i+1).zfill(2)+'.nc'

		if year_i == 300 and month_i == 0:
			filename = directory_data+'t.t0.1_42l_nccs01.'+str(301).zfill(4)+str(month_i+1).zfill(2)+'.nc'
			
		print(filename)
		lon, lat, depth, layer, depth_top, depth_u, grid_x, v_vel_month = ReadinData(filename, lat_index)

		v_vel[month_i]	= v_vel_month

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	v_vel	= np.sum(v_vel * month_days_all, axis = 0)

	#------------------------------------------------------------------------------

	#Determine the meridional transport
	transport	= v_vel * layer_field_regions * grid_x

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[year_i - int(np.min(time))]	= np.sum(transport) / 1000000.0
	
#-----------------------------------------------------------------------------------------

plot(time_year, transport_all, '-r')
show()


print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'w')

fh.createDimension('time', len(time_year))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)

fh.variables['Transport'].long_name 	= 'Volume transport'

fh.variables['time'].units 		= 'Year'
fh.variables['Transport'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time_year
fh.variables['Transport'][:]    	= transport_all

fh.close()
