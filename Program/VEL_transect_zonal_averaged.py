#Program determines zonal and annual averaged zonal velocity

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
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0

	#Select region
	lon_min_index	= (fabs(lon - -60)).argmin()
	lon_max_index	= (fabs(lon - 0)).argmin() + 1
	lat_min_index	= (fabs(lat - -75)).argmin()
	lat_max_index	= (fabs(lat - 10)).argmin() + 1	
	
	lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
		
	fh = netcdf.Dataset(filename, 'r')

	u_vel 		= fh.variables['UVEL'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index] / 100.0	#Zonal velocity (m/s)
	v_vel 		= fh.variables['VVEL'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index] / 100.0	#Meridional velocity (m/s)
	
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
			u_vel_depth		= ma.masked_where((u_vel_depth <= 0) & (v_vel_depth <= 0), u_vel_depth) #mask cells where both uvel and vvel are zero (then it is some kind of bug in the HR pop)
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
		volume			= ma.masked_where(volume <= 0, volume)
		volume_norm		= ma.masked_all(np.shape(volume))
		
		for depth_i in range(len(depth)):
			for lat_i in range(len(lat)):	
				#Normalise the field for each latitude and depth layer
				volume_norm[depth_i, lat_i]	= volume[depth_i, lat_i] / np.sum(volume[depth_i, lat_i])
				
	#Take the volume-averaged zonal mean for u_vel
	u_vel		= np.sum(u_vel * volume_norm, axis = 2)
	
	print(np.min(u_vel))
	print(np.max(u_vel))
	
	#CS	= contourf(lat, depth, u_vel, levels = np.arange(-0.5, 0.51, 0.1), extend = 'both', cmap = 'bwr')
	#colorbar(CS)
	#show()
	#sys.exit()
		
	return lat, depth, volume_norm, u_vel
	
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#Second SOM cycle (model year 45-91), last SOM cycle (495-571)
year_start	= 495
year_end	= 571

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

#note that january 300 does not exist and that february 300 does not have PD data. Therefore, we take january 301 to january 351
files	= files[year_start*12-1:year_end*12]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------
lat, depth, volume_norm, u_vel		= ReadinData(files[0])
time_year				= ma.masked_all(year_end-year_start+1)
u_vel_all				= ma.masked_all((len(time_year), len(depth), len(lat)))

u_vel_transect				= ma.masked_all((len(depth), len(lat)))

print(np.shape(u_vel))
print(np.shape(depth))
print(np.shape(lat))

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	u_vel_year 		= ma.masked_all((12, len(depth), len(lat)))
		
	for month_i in range(12):

		#Get the monthly files
		filename 	= files[year_i*12 + month_i]

		print(filename)
		lat, depth, volume_norm, u_vel_year[month_i] = ReadinData(filename)

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lat)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	u_vel_all[year_i]	= np.sum(u_vel_year * month_days_all, axis = 0)
	
#Annual mean
u_vel_transect = np.nanmean(u_vel_all, axis = 0)
	
print(np.shape(u_vel_all))

contourf(lat, depth, u_vel_transect)
show()	

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/UVEL_year_'+str(year_start)+'-'+str(year_end)+'_zonal_averaged_60W_0W_transect_SO.nc', 'w')

fh.createDimension('lat', len(lat))
fh.createDimension('depth', len(depth))

fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('depth', float, ('depth'), zlib=True)
fh.createVariable('U_VEL', float, ('depth', 'lat'), zlib=True)

fh.variables['lat'].long_name 		= 'Array of u-latitudes'
fh.variables['depth'].long_name 	= 'Depth'
fh.variables['U_VEL'].long_name 	= 'Zonal velocity'

fh.variables['depth'].units 		= 'm'
fh.variables['lat'].units 		= 'Degrees N'
fh.variables['U_VEL'].units 		= 'm/s'

#Writing data to correct variable
fh.variables['depth'][:] 		= depth
fh.variables['lat'][:] 			= lat
fh.variables['U_VEL'][:] 		= u_vel_transect

fh.close()	
	
	
