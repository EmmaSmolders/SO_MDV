#Non-zonality parameter 

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
	
	#print(lat)
	
	#plt.figure()
	#plt.plot(lat)
	#plt.show()
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0
	
	#plt.figure()
	#plt.plot(lon)
	#plt.show()

	#print(lon)
	#print(np.shape(lon))
	#print(np.min(lon))

	#Select region 
	#lon_min_index	= (fabs(lon - 0)).argmin()
	#lon_max_index	= (fabs(lon - 40)).argmin() + 1
	#lat_min_index	= (fabs(lat - -90)).argmin()
	#lat_max_index	= (fabs(lat - -45)).argmin() + 1
	#lon_min_index	= (fabs(lon - -50)).argmin()
	#lon_max_index	= (fabs(lon - 0)).argmin() + 1
	lat_min_index	= (fabs(lat - -90)).argmin()
	lat_max_index	= (fabs(lat - -30)).argmin() + 1
	depth_min_index = (fabs(depth - 5)).argmin()
	depth_max_index = (fabs(depth - 318)).argmin() + 1
	
	#Test region
	#lon_min_index	= (fabs(lon - 0)).argmin()
	#lon_max_index	= (fabs(lon - 1)).argmin() + 1
	#lat_min_index	= (fabs(lat - -47)).argmin()
	#lat_max_index	= (fabs(lat - -45)).argmin() + 1
	#depth_min_index = (fabs(depth - 5)).argmin()
	#depth_max_index = (fabs(depth - 20)).argmin() + 1
	
	#print(lon_min_index)
	#print(lon_max_index)
	#print(lat_min_index)
	#print(lat_max_index)
	
	#sys.exit()
	
	#print(depth_min_index)
	#print(depth_max_index)	
	
	#lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, :]
	area		= area[lat_min_index:lat_max_index, :]
	depth		= depth[depth_min_index:depth_max_index]
	layer		= layer[depth_min_index:depth_max_index]
	
	print(lon)
	print(lat)
	print(depth)
	
	#sys.exit()
		
	fh = netcdf.Dataset(filename, 'r')

	u_vel 		= fh.variables['UVEL'][depth_min_index:depth_max_index, lat_min_index:lat_max_index, :] / 100.0	#Zonal velocity (m/s)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_min_index:lat_max_index, :] / 100.0	#Meridional velocity (m/s)
	
	#print(np.shape(layer))
	
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
		
		#print(np.shape(volume))
		#print(np.shape(u_vel))
		
		#Volume integrated velocities
		U_vel 			= np.sum(u_vel*volume)
		V_vel			= np.sum(v_vel*volume)
		
		#print(U_vel)
		
	return lat, lon, depth, volume, u_vel, v_vel
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#Second SOM cycle (model year 63-114), second (324-378) or last SOM cycle (500-600)
year_start	= 324
year_end	= 378

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

window = 5

#-----------------------------------------------------------------------------------------

lat, lon, depth, volume, u_vel, v_vel		= ReadinData(files[0])

print(np.shape(volume))

#sys.exit()

print(depth)

time_year					= ma.masked_all(year_end-year_start+1)
#VVUU_vel_all					= ma.masked_all((len(time_year)))
uu_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
vv_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
#vvuu_vel_all					= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
V2U2_int  					= ma.masked_all((len(time_year) - window + 1))

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	for month_i in range(12):

		#Get the monthly files
		filename 	= files[year_i*12 + month_i]

		#print(filename)
		lat, lon, depth, volume, u_vel_year, v_vel_year = ReadinData(filename)
		
		# Compute the rolling index in uu_vel_all for the current month
		pdw_index = (year_i % window) * 12 + month_i
		#print(f"Year: {year_i}, Month: {month_i}, PDW_all index: {pdw_index}")
		
		#Save monthly u*u in window
		uu_vel_all[pdw_index,:,:,:] = u_vel_year*u_vel_year
		vv_vel_all[pdw_index,:,:,:] = v_vel_year*v_vel_year
		
	#Start calculating energetics when year_i has at least 5 years
	if year_i >= window - 1:

		#5-year window available
		V2 = np.mean(vv_vel_all, axis=0)
		U2 = np.mean(uu_vel_all, axis=0)
		V2U2 = V2/U2

		#Volume integratedv
		V2U2_int[year_i - window + 1] = np.sum(V2U2*volume)
	
	#-----------------------------------------------------------------------------------------
	#-----------------------------------------------------------------------------------------

	print('Data is written to file')
	fh = netcdf.Dataset(directory+'Ocean/Nonzonality_year_'+str(year_start)+'-'+str(year_end)+'_volume_integrated_SO30AREA_window_'+str(window)+'.nc', 'w')

	fh.createDimension('time_int', len(time_year) - window + 1)
	fh.createDimension('time', len(time_year))
	fh.createDimension('lat', len(lat))
	fh.createDimension('lon', len(lon))

	fh.createVariable('time', float, ('time'), zlib=True)
	fh.createVariable('time_int', float, ('time_int'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	#fh.createVariable('V2_U2', float, ('time'), zlib=True)
	fh.createVariable('V2_U2_int', float, ('time_int'), zlib=True)

	fh.variables['time'].long_name 		= 'Model year'
	fh.variables['time_int'].long_name 	= 'Starting year of window'
	fh.variables['lat'].long_name		= 'Latitudes for surface integration'
	fh.variables['lon'].long_name		= 'Longitudes for surface integration'
	#fh.variables['V2_U2'].long_name 	= 'Volume integrated nonzonality parameter (V^2/U^2). First integrated, then division'
	fh.variables['V2_U2_int'].long_name 	= 'Volume integrated nonzonality parameter (V^2/U^2)'

	fh.variables['time'].units 		= 'model year'
	fh.variables['time_int'].units 		= 'model year'
	fh.variables['lat'].units 		= 'Degrees N'
	fh.variables['lon'].units 		= 'Degrees E'
	#fh.variables['V2_U2'].units 		= ''
	fh.variables['V2_U2_int'].units 	= ''

	#Writing data to correct variable
	fh.variables['time'][:] 	= time_year
	fh.variables['time_int'][:] 	= time_year[:len(time_year) - window + 1]
	fh.variables['lat'][:] 		= lat
	fh.variables['lon'][:] 		= lon
	#fh.variables['V2_U2'][:] 	= VVUU_vel_all
	fh.variables['V2_U2_int'][:] 	= V2U2_int

	fh.close()









