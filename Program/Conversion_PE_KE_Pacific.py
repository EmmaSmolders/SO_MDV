#Program determines conversion PE to KE of the SOM (following Juling's code)

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory_data		= '/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/tavg/'
directory_data_wvel	= '/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/wvel/'
directory		= '/home/smolders/HR_POP/Data/'

def ReadinData(filename, layer_avail = False, volume_norm = False):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#First get the u-grid
	lon 		= fh.variables['TLONG'][400]				#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['TLAT'][:, 780]				#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	dx		= fh.variables['DXT'][:]  / 100.0			#Zonal grid cell length (m)
	depth  	 	= fh.variables['z_t'][:]  / 100.0			#Depth (m)
	depth_grid 	= fh.variables['HT'][:] / 100.0				#Depth at t-grid (m)
	layer		= fh.variables['dz'][:] / 100.0				#Layer thickness (m)
	depth_top	= fh.variables['z_w_top'][:] / 100.0			#Top of grid cell
	area		= fh.variables['TAREA'][:] / 10000			#TAREA (m^2)

	fh.close()
	
	#lon must be from 0 to 360
	lon[lon<0]	= lon[lon<0] + 360.0
	
	#Negative longitudes are also allowed
	lon_min 	= 150
	lon_max	        = 250
	lat_min		= -70
	lat_max		= -45
	
	#First get the latitude extent
	lat_min_index	= (fabs(lat - lat_min)).argmin()
	lat_max_index	= (fabs(lat - lat_max)).argmin() + 1	
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid 	= depth_grid[lat_min_index:lat_max_index,:]
	area		= area[lat_min_index:lat_max_index,:]
	
	#sys.exit()
		
	fh = netcdf.Dataset(filename, 'r')

	PD 		= fh.variables['PD'][:, lat_min_index:lat_max_index, :]*1000 - 1026	#Potential density (kg/m^3) - rho0
	
	fh.close()
	
	#print(np.where(PD < 0))
	#print(np.nanmean(PD))
	
	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		PD[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], PD[depth_i])
		
	#Convert the field to 0E to 360E
	lon_2, area	= ConverterField2D(lon, area)
	lon_2, PD	= ConverterField3D(lon, PD)
	lon, depth_grid	= ConverterField2D(lon, depth_grid)
	
	print(lon)

	if lon_min < 0 or lon_max > 360:
		#If needed, extend the field to add negative longitudes and beyond 360E
		lon_2, PD		= PeriodicBoundaries3D(lon, lat, PD, 1800)
		lon_2, depth_grid	= PeriodicBoundaries2D(lon, lat, depth_grid, 1800)
		lon, area		= PeriodicBoundaries2D(lon, lat, area, 1800)	
		
	print(lon)
	
	#Get the longitude extent
	lon_min_index	= (fabs(lon - lon_min)).argmin()
	lon_max_index	= (fabs(lon - lon_max)).argmin() + 1	
	lon		= lon[lon_min_index:lon_max_index]
	area		= area[:, lon_min_index:lon_max_index]
	PD		= PD[:,:, lon_min_index:lon_max_index]
	depth_grid	= depth_grid[:, lon_min_index:lon_max_index]
	
	print(lon)
	print(lat)

	#------------------------------------------------------------------------------
	#Make volume on T-grid
	if layer_avail	== False:	
		#Determine the depth per grid cell (parcel bottom cells)
		layer_field	= ma.masked_all((len(depth), len(lat), len(lon)))

		for depth_i in range(len(layer)):
			#print(depth_i)

			#Mask all elements which are land and fill the layer field with the depth layer for each layer
			PD_depth		= PD[depth_i]
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
		volume			= ma.masked_where(volume <= 0, volume) #volume of gridcells
		
	return lat, lon, depth, area, layer, volume, PD
	
def ConverterField2D(lon, field):
	"""Shifts field"""
	lon_new		= ma.masked_all(shape(lon))
	field_new	= ma.masked_all(shape(field))

	#Get the corresponding index
	index		= (fabs(lon - 0)).argmin()

	#Start filling at -180
	lon_new[:len(lon[index:])] 		= lon[index:]
	field_new[:, :len(lon[index:])]		= field[:, index:]

	#Fill the remaining part
	lon_new[len(lon[index:]):] 		= lon[:index]
	field_new[:, len(lon[index:]):]		= field[:, :index]

	return lon_new, field_new

def PeriodicBoundaries2D(lon, lat, field, lon_grids = 1):
        """Add periodic zonal boundaries for 2D field"""

        #Empty field with additional zonal boundaries
        lon_2                   = np.zeros(len(lon) + lon_grids * 2)
        field_2                 = ma.masked_all((len(lat), len(lon_2)))

        #Get the left boundary, which is the right boundary of the original field
        lon_2[:lon_grids]       = lon[-lon_grids:] - 360.0
        field_2[:, :lon_grids]  = field[:, -lon_grids:]

        #Same for the right boundary
        lon_2[-lon_grids:]      = lon[:lon_grids] + 360.0
        field_2[:, -lon_grids:] = field[:, :lon_grids]

        #And the complete field
        lon_2[lon_grids:-lon_grids]             = lon
        field_2[:, lon_grids:-lon_grids]        = field

        return lon_2, field_2

def ConverterField3D(lon, field):
	"""Shifts field, to -180E to 180E"""
	lon_new		= ma.masked_all(shape(lon))
	field_new	= ma.masked_all(shape(field))

	#Get the corresponding index
	index		= (fabs(lon - 0)).argmin()

	#Start filling at -180E
	lon_new[:len(lon[index:])] 	 	= lon[index:]
	field_new[:, :, :len(lon[index:])]	= field[:, :, index:]

	#Fill the remaining part to 360E
	lon_new[len(lon[index:]):]	 	= lon[:index]
	field_new[:, :, len(lon[index:]):]	= field[:, :, :index]

	#lon_new[lon_new >= 179.9]	= lon_new[lon_new >= 179.9] - 360.0

	return lon_new, field_new
	
def PeriodicBoundaries3D(depth, lon, lat, field, lon_grids = 1):
        """Add periodic zonal boundaries for 3D field"""

        #Empty field with additional zonal boundaries
        lon_2                   = np.zeros(len(lon) + lon_grids * 2)
        field_2                 = ma.masked_all((len(depth), len(lat), len(lon_2)))
        
        print(np.shape(field_2))
        print(len(lat))
        print(len(lon_2))

        #Get the left boundary, which is the right boundary of the original field
        lon_2[:lon_grids]       = lon[-lon_grids:] - 360.0
        field_2[:, :, :lon_grids]  = field[:, :, -lon_grids:]

        #Same for the right boundary
        lon_2[-lon_grids:]      = lon[:lon_grids] + 360.0
        field_2[:, :, -lon_grids:] = field[:, :, :lon_grids]

        #And the complete field
        lon_2[lon_grids:-lon_grids]             = lon
        field_2[:, :, lon_grids:-lon_grids]     = field

        return lon_2, field_2
	
def ReadinDataWVEL(filename, layer_avail = False, volume_norm = False):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#First get the u-grid
	lon 		= fh.variables['TLONG'][400]				#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['TLAT'][:, 780]				#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	dx		= fh.variables['DXT'][:]  / 100.0			#Zonal grid cell length (m)
	depth  	 	= fh.variables['z_t'][:]  / 100.0			#Depth (m)
	depth_grid 	= fh.variables['HT'][:] / 100.0				#Depth at U-grid (m)
	layer		= fh.variables['dz'][:] / 100.0				#Layer thickness (m)
	depth_top	= fh.variables['z_w_top'][:] / 100.0			#Top of grid cell
	area		= fh.variables['TAREA'][:] / 10000			#UAREA (m^2)

	fh.close()
	
	#lon must be from 0 to 360
	lon[lon<0]	= lon[lon<0] + 360.0
	
	#Negative longitudes are also allowed
	lon_min 	= 150
	lon_max	        = 250
	lat_min		= -70
	lat_max		= -45
	
	#First get the latitude extent
	lat_min_index	= (fabs(lat - lat_min)).argmin()
	lat_max_index	= (fabs(lat - lat_max)).argmin() + 1	
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid 	= depth_grid[lat_min_index:lat_max_index,:]
	area		= area[lat_min_index:lat_max_index,:]
	
	#sys.exit()
		
	fh = netcdf.Dataset(filename, 'r')

	w_vel 		= fh.variables['WVEL'][:, lat_min_index:lat_max_index, :]/ 100.0	#Vertical velocity (m/s)
	
	fh.close()
	
	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		w_vel[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], w_vel[depth_i])
		
	#Convert the field to 0E to 360E
	lon, w_vel	= ConverterField3D(lon, w_vel)
	
	print(lon)

	if lon_min < 0 or lon_max > 360:
		#If needed, extend the field to add negative longitudes and beyond 360E
		lon, w_vel		= PeriodicBoundaries3D(lon, lat, w_vel, 1800)
		
	print(lon)
	
	#Get the longitude extent
	lon_min_index	= (fabs(lon - lon_min)).argmin()
	lon_max_index	= (fabs(lon - lon_max)).argmin() + 1	
	lon		= lon[lon_min_index:lon_max_index]
	w_vel		= w_vel[:,:, lon_min_index:lon_max_index]
	
	print(lon)
	print(lat)
		
	return w_vel

def wtt2ttt(W_FIELD, DZT):
	"""
    	Interpolates a field from the WTT-grid to the TTT-grid.

    	Parameters:
        	W_FIELD (numpy.ndarray): Input field on the WTT-grid (shape: [depth, lat, lon]).
        	DZT (numpy.ndarray): Layer thickness on the TTT-grid (shape: [depth, lat, lon]).

    	Returns:
        	numpy.ndarray: Interpolated field on the TTT-grid (shape: [depth, lat, lon]).
    	"""
	
	#Get dimensions
	km, jmt, imt = W_FIELD.shape
	print(km)
	print(jmt)
	print(imt)
	
	#Initialize T_FIELD with zeros
	T_FIELD = ma.masked_all((km, jmt, imt))
	
	#Interpolate for all levels except bottom level
	for k in range(km -1):
		if DZT[k] != 0:
			T_FIELD[k, :, :] = 0.5 * (W_FIELD[k, :, :] + W_FIELD[k + 1, : , :])
			
	#Bottom level
	if DZT[km-1] != 0:
		T_FIELD[km-1, :, :] = 0.5 * W_FIELD[k-1, :, :]
		
	return T_FIELD
			
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#First SOM cycle (model year 63-114), second (324-378) or last SOM cycle (500-600)
year_start	= 63
year_end	= 600

window = 5
window_counter = 0

#-----------------------------------------------------------------------------------------

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

#note that january 300 is 300-02 etc. and that february 300 does not have PD data. 
files	= files[(year_start-1)*12:year_end*12]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------

files_wvel = []

files_wvel_all	= glob.glob(directory_data_wvel+'WVEL_t.t0.1_42l_nccs01.*.nc')
files_wvel_all.sort()

print(len(files_wvel_all[0]))

for file_i in range(len(files_wvel_all)):

	if len(files_wvel_all[file_i]) == 121:
		files_wvel.append(files_wvel_all[file_i])

#february 300 does not have PD data. 
files_wvel	= files_wvel[(year_start-1)*12:year_end*12]

print(files_wvel[0])
print(files_wvel[-1])

#-----------------------------------------------------------------------------------------
lat, lon, depth, area, layer, volume, PD	= ReadinData(files[0])
w_vel						= ReadinDataWVEL(files_wvel[0])

time_year			= ma.masked_all(year_end-year_start+1)
PD_all				= ma.masked_all((window, len(depth), len(lat), len(lon)))
w_vel_all			= ma.masked_all((window, len(depth), len(lat), len(lon)))
PDW_all				= ma.masked_all((window*12, len(depth), len(lat), len(lon)))

#Volume integrated energetics
C_mean_int  	= ma.masked_all((len(time_year) - window + 1))
C_eddy_int  	= ma.masked_all((len(time_year) - window + 1))
C_total_int  	= ma.masked_all((len(time_year) - window + 1))
	
#Gravitational constant
g = 9.81
rho0 = 1026

window_counter = 0

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start
	#print(time_year[year_i])
	
	PD_year 		= ma.masked_all((12, len(depth), len(lat), len(lon)))
	w_vel_year 		= ma.masked_all((12, len(depth), len(lat), len(lon)))

	for month_i in range(12):
		
		#print(year_i+300+year_start-1)
		
		#Get the monthly files 
		filename 	= files[year_i*12 + month_i]
		#print(filename)
		
		filename_w_vel 	= files_wvel[year_i*12 + month_i]
		#print(filename_w_vel)
		
		lat, lon, depth, area, layer, volume, PD_year[month_i]	 = ReadinData(filename)
		w_vel_year[month_i] 		 = ReadinDataWVEL(filename_w_vel)
		
		# Compute the rolling index in PDW_all for the current month
		pdw_index = (year_i % window) * 12 + month_i
        	
		#print(f"Year: {year_i}, Month: {month_i}, PDW_all index: {pdw_index}")
		
		for k in range(len(depth)):
			if k == 0:
				PDW_all[pdw_index, k, :,:] = 0.5*w_vel_year[month_i,k+1,:,:] * PD_year[month_i,k,:,:]
			if k == len(depth) - 1:
				PDW_all[pdw_index, k, :,:] = 0.5*w_vel_year[month_i,k,:,:] * PD_year[month_i,k,:,:]
			else:
				PDW_all[pdw_index, k, :,:] = 0.5*(w_vel_year[month_i,k,:,:] + w_vel_year[month_i,k+1,:,:]) * PD_year[month_i,k,:,:]

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lat), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	# Determine the time mean over the months of choice in the window length (sliding window mean)
	w_vel_all[year_i % window] 	= np.sum(w_vel_year * month_days_all, axis=0)
	PD_all[year_i % window] 	= np.sum(PD_year * month_days_all, axis=0)
	
	#Start calculating energetics when year_i has at least 5 years
	if year_i >= window - 1:
	
		#Compute mean values in window
		w_vel_mean 	= np.mean(w_vel_all, axis = 0)
		PD_mean		= np.mean(PD_all, axis = 0)
		PDW_mean	= np.mean(PDW_all, axis = 0)
		
		print(np.mean(w_vel_mean))
		print(np.mean(PD_mean))
		
		TTT_PDW = wtt2ttt(PDW_mean, layer)
		cPKt = -g * TTT_PDW
		
		cPKm = np.zeros_like(w_vel_mean)
		cPKe = np.zeros_like(w_vel_mean)
		
		# 5-year window available
		for k in range(len(depth)):
    			if k == 0:  # First depth level
        			cPKm[k, :, :] = w_vel_mean[k+1,:,:] * (PD_mean[k, :, :] + PD_mean[k+1, :, :])
    			elif k == len(depth) - 1:  # Last depth level
        			cPKm[k, :, :] = w_vel_mean[k,:,:] * (PD_mean[k, :, :] + PD_mean[k-1, :, :])
    			else:  # Intermediate depth levels
        			cPKm[k, :, :] = w_vel_mean[k+1,:,:] * (PD_mean[k, :, :] + PD_mean[k+1, :, :]) + w_vel_mean[k,:,:] * (PD_mean[k, :, :] + PD_mean[k-1, :, :])
    
    			cPKm[k, :, :] = -g * 0.25 * cPKm[k, :, :]
    			cPKe[k, :, :] = cPKt[k, :, :] - cPKm[k, :, :]
    			
		#print(np.mean(cPKm))
		#print(np.mean(cPKt))
		#print(np.mean(cPKe))

		# Volume integrated
		C_mean_int[year_i - window + 1] 	= np.sum(cPKm * volume)
		C_total_int[year_i - window + 1] 	= np.sum(cPKt * volume)
		C_eddy_int[year_i - window + 1] 	= np.sum(cPKe * volume)
		
		#print(C_mean_int[year_i - window + 1])
		#print(C_eddy_int[year_i - window + 1])
		#sys.exit()

	#-----------------------------------------------------------------------------------------
	#-----------------------------------------------------------------------------------------

	print('Data is written to file')
	fh = netcdf.Dataset(directory+'Ocean/Conversion_PE_KE_year_'+str(year_start)+'-'+str(year_end)+'_'+str(window)+'_volume_integrated_Pacific_testjuling.nc', 'w')

	fh.createDimension('time', len(time_year) - window + 1)
	fh.createDimension('lat', len(lat))
	fh.createDimension('lon', len(lon))

	fh.createVariable('time', float, ('time'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('C_mean', float, ('time'), zlib=True)
	fh.createVariable('C_total', float, ('time'), zlib=True)
	fh.createVariable('C_eddy', float, ('time'), zlib=True)

	fh.variables['time'].long_name 		= 'Starting year of window'
	fh.variables['lat'].long_name		= 'Latitudes for surface integration'
	fh.variables['lon'].long_name		= 'Longitudes for surface integration'
	fh.variables['C_total'].long_name 	= 'Volume integrated total conversion energy'
	fh.variables['C_mean'].long_name 	= 'Volume integrated mean conversion energy'
	fh.variables['C_eddy'].long_name 	= 'Volume integrated eddy conversion energy'

	fh.variables['time'].units 	= 'year'
	fh.variables['lat'].units 	= 'Degrees N'
	fh.variables['lon'].units 	= 'Degrees E'
	fh.variables['C_mean'].units 	= 'W'
	fh.variables['C_total'].units 	= 'W'
	fh.variables['C_eddy'].units 	= 'W'

	#Writing data to correct variable
	fh.variables['time'][:] 	= time_year[:len(time_year) - window + 1]
	fh.variables['lat'][:] 		= lat
	fh.variables['lon'][:] 		= lon
	fh.variables['C_total'][:] 	= C_total_int
	fh.variables['C_mean'][:] 	= C_mean_int
	fh.variables['C_eddy'][:] 	= C_eddy_int

	fh.close()	
	
	
	
	
	
	
	
