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
	lon 		= fh.variables['TLONG'][400]				#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['TLAT'][:, 780]				#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	dx		= fh.variables['DXT'][:]  / 100.0			#Zonal grid cell length (m)
	depth  	 	= fh.variables['z_t'][:]  / 100.0			#Depth (m)
	depth_grid 	= fh.variables['HT'][:] / 100.0				#Depth at t-grid (m)
	layer		= fh.variables['dz'][:] / 100.0				#Layer thickness (m)
	depth_top	= fh.variables['z_w_top'][:] / 100.0			#Top of grid cell
	area		= fh.variables['TAREA'][:] / 10000			#TAREA (m^2)

	fh.close()
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0

	#Select region (SO30 or WGKP region)
	lon_min_index	= (fabs(lon - -35)).argmin()
	lon_max_index	= (fabs(lon - 80)).argmin() + 1
	lat_min_index	= (fabs(lat - -90)).argmin()
	lat_max_index	= (fabs(lat - -50)).argmin() + 1
	
	lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	
	#Select region (SO30 or WGKP region)
	#lat_min_index	= (fabs(lat - -90)).argmin()
	#lat_max_index	= (fabs(lat - -30)).argmin() + 1
	
	#lon		= lon[lon_min_index:lon_max_index]
	#lat		= lat[lat_min_index:lat_max_index]
	#depth_grid	= depth_grid[lat_min_index:lat_max_index, :]
	#area		= area[lat_min_index:lat_max_index, :]
	
	#print(lon)
	#print(lat)
	
	#sys.exit()
		
	fh = netcdf.Dataset(filename, 'r')

	PD 		= fh.variables['PD'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index]*1000 	#Potential density (kg/m^3)
	
	fh.close()
	
	print(np.where(PD < 0))
	
	for depth_i in range(len(depth)):
		#Mask all the field at the topography
		PD[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], PD[depth_i])

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
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

# Constants
g = 9.81  # Gravitational constant
window = 5 # Moving average window

#First SOM cycle (model year 63-114), second (324-378) or last SOM cycle (500-600)
year_start	= 324
year_end	= 378

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

#note that january 300 is written in feb 300. (so 300-01 == 300-02 etc.)
files	= files[(year_start-1)*12:year_end*12]

print(files[0])
print(files[-1])

# Read initial data
lat, lon, depth, area, layer, volume, PD = ReadinData(files[0])

# Initialize arrays
time_year = np.arange(year_start, year_end+1)

# Precompute month weights
month_days = np.array([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
month_weights = month_days / np.sum(month_days)

# Process each year
for year_i, year in enumerate(time_year):
    PD_year = ma.masked_all((12, len(depth), len(lat), len(lon)))
    for month_i in range(12):
        filename = files[year_i * 12 + month_i]
        _, _, _, _, _, _, PD_year[month_i] = ReadinData(filename)

    # Compute weighted mean for the year
    PD_all = np.sum(PD_year * month_weights[:, np.newaxis, np.newaxis, np.newaxis], axis=0)

    if year_i == 0:
        PD_mean = np.copy(PD_all) / len(time_year) 

    else:
        PD_mean += np.copy(PD_all) / len(time_year) 
            
# Compute reference state
#PD_mean = np.mean(PD_all, axis=0)

#del PD_all

rho_ref = np.sum(PD_mean * area[np.newaxis, :, :], axis=(1, 2)) / np.sum(area)

del PD_mean

# Compute vertical gradient of rho_ref
n0 = np.zeros(len(depth))
n0[1:-1] = (rho_ref[:-2] - rho_ref[2:]) / (2.0 * layer[1:-1])
n0[0] = (rho_ref[0] - rho_ref[1]) / layer[0]
n0[-1] = (rho_ref[-2] - rho_ref[-1]) / layer[-1]

del rho_ref

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/n0_year_'+str(year_start)+'-'+str(year_end)+'_SO30.nc', 'w')

fh.createDimension('depth', len(depth))
fh.createDimension('lat', len(lat))
fh.createDimension('lon', len(lon))

fh.createVariable('depth', float, ('depth'), zlib=True)
fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('lon', float, ('lon'), zlib=True)
fh.createVariable('n0', float, ('depth'), zlib=True)

fh.variables['depth'].long_name 	= 'Depth'
fh.variables['lat'].long_name		= 'Latitudes'
fh.variables['lon'].long_name		= 'Longitudes'
fh.variables['n0'].long_name 		= 'Inverse vertical density gradient'

fh.variables['depth'].units 	= 'm'
fh.variables['lat'].units 	= 'Degrees N'
fh.variables['lon'].units 	= 'Degrees E'
fh.variables['n0'].units 	= 'kg/m^4'

#Writing data to correct variable
fh.variables['depth'][:] 	= depth
fh.variables['lat'][:] 		= lat
fh.variables['lon'][:] 		= lon
fh.variables['n0'][:] 		= n0

fh.close()

#Read in mean PD
#fh = netcdf.Dataset(directory+'Ocean/PD_mean_year_'+str(year_start)+'-'+str(year_end)+'_SO30.nc','r')

#PD_mean = fh.variables['PD_mean'][:] #time mean PD (shape [depth, lat, lon])

#fh.close()

#rho_ref = np.sum(PD_mean * area[np.newaxis, :, :], axis=(1, 2)) / np.sum(area)

#del PD_mean

# Compute vertical gradient of rho_ref
#n0 = np.zeros(len(depth))
#n0[1:-1] = (rho_ref[:-2] - rho_ref[2:]) / (2.0 * layer[1:-1])
#n0[0] = (rho_ref[0] - rho_ref[1]) / layer[0]
#n0[-1] = (rho_ref[-2] - rho_ref[-1]) / layer[-1]

#del rho_ref

# Initialize arrays
time_year 	= np.arange(year_start, year_end+1)
PDPD_all 	= ma.masked_all((window*12, len(depth), len(lat), len(lon)))
APE_int 	= ma.masked_all(len(time_year) - window + 1)

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start
	#print(time_year[year_i])

	for month_i in range(12):
		
		#print(year_i+300+year_start-1)
		
		#Get the monthly files 
		filename 	= files[year_i*12 + month_i]
		print(filename)
		
		lat, lon, depth, area, layer, volume, PD_month	 = ReadinData(filename)
		
		# Compute the rolling index in PDW_all for the current month
		pdw_index = (year_i % window) * 12 + month_i
		#print(f"Year: {year_i}, Month: {month_i}, PDW_all index: {pdw_index}")
		
		#Save monthly PD*PD in window
		PDPD_all[pdw_index, :,:,:] = PD_month**2
        	
	#Start calculating energetics when year_i has at least 5 years
	if year_i >= window - 1:
	
		#Compute mean values in window
		PD2_mean	= np.mean(PDPD_all, axis = 0)
		
		APE_year = -0.5 * g * (1 / n0[:, np.newaxis, np.newaxis]) * PD2_mean
		
		#Volume integrated
		APE_int[year_i - window + 1] 	= np.sum(APE_year*volume)
		
	#-----------------------------------------------------------------------------------------
	#-----------------------------------------------------------------------------------------

	print('Data is written to file')
	fh = netcdf.Dataset(directory+'Ocean/APE_year_'+str(year_start)+'-'+str(year_end)+'_window_'+str(window)+'_volume_integrated_WGKP.nc', 'w')

	fh.createDimension('time', len(time_year) - window + 1)
	fh.createDimension('lat', len(lat))
	fh.createDimension('lon', len(lon))

	fh.createVariable('time', float, ('time'), zlib=True)
	fh.createVariable('lat', float, ('lat'), zlib=True)
	fh.createVariable('lon', float, ('lon'), zlib=True)
	fh.createVariable('APE', float, ('time'), zlib=True)

	fh.variables['time'].long_name 	= 'Starting year of window'
	fh.variables['lat'].long_name	= 'Latitudes for surface integration'
	fh.variables['lon'].long_name	= 'Longitudes for surface integration'
	fh.variables['APE'].long_name 	= 'Volume integrated available potential energy'

	fh.variables['time'].units 	= 'year'
	fh.variables['lat'].units 	= 'Degrees N'
	fh.variables['lon'].units 	= 'Degrees E'
	fh.variables['APE'].units 	= 'J'

	#Writing data to correct variable
	fh.variables['time'][:] 	= time_year[:len(time_year) - window + 1]
	fh.variables['lat'][:] 		= lat
	fh.variables['lon'][:] 		= lon
	fh.variables['APE'][:] 		= APE_int

	fh.close()	
		
