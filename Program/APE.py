from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import gsw

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

	#Select region (WGKP region)
	#lon_min_index	= (fabs(lon - -35)).argmin()
	#lon_max_index	= (fabs(lon - 80)).argmin() + 1
	#lat_min_index	= (fabs(lat - -90)).argmin()
	#lat_max_index	= (fabs(lat - -50)).argmin() + 1
	
	#lon		= lon[lon_min_index:lon_max_index]
	#lat		= lat[lat_min_index:lat_max_index]
	#depth_grid	= depth_grid[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	#area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	
	#Select region (SO30 region)
	lat_min_index	= (fabs(lat - -90)).argmin()
	lat_max_index	= (fabs(lat - -30)).argmin() + 1
	
	#lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_grid	= depth_grid[lat_min_index:lat_max_index, :]
	area		= area[lat_min_index:lat_max_index, :]
	
	#print(lon)
	#print(lat)
	
	#sys.exit()
		
	fh = netcdf.Dataset(filename, 'r')

	#PD 		= fh.variables['PD'][:, lat_min_index:lat_max_index, lon_min_index:lon_max_index]*1000 	#Potential density (kg/m^3)
	#PD 		= fh.variables['PD'][:, lat_min_index:lat_max_index, :]*1000 	#Potential density (kg/m^3)
	temp 		= fh.variables['TEMP'][:, lat_min_index:lat_max_index, :] 	#Potential temperature (degC)
	salt 		= fh.variables['SALT'][:, lat_min_index:lat_max_index, :]*1000 	#Salinity (g/kg)
	
	fh.close()
	
	#Mask the area
	temp		= ma.masked_where(salt <= 0, temp)	
	salt		= ma.masked_where(salt <= 0, salt)
	
	PD 		= ma.masked_all((len(depth), len(lat), len(lon)))
	
	for depth_i in range(len(depth)):
		#print(depth_i)
		#First determine the conservative temperature from the potential temperature
		temp_CT		= gsw.CT_from_pt(salt[depth_i], temp[depth_i])
		
		#Get the potential density
		PD[depth_i]		= gsw.sigma0(salt[depth_i], temp_CT)+1000.0
	
	#print(np.where(PD < 0))
	
	#for depth_i in range(len(depth)):
		#Mask all the field at the topography
	#	PD[depth_i]	= ma.masked_where(depth_grid <= depth_top[depth_i], PD[depth_i])

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

#First SOM cycle (model year 63-114), second (410-480) or last SOM cycle (500-600)
year_start	= 350
year_end	= 600

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

#-----------------------------------------------------------------------------------------

# Read initial data
lat, lon, depth, area, layer, volume, PD = ReadinData(files[0])

# Read in reference state and take time average
fh = netcdf.Dataset(directory + 'Ocean/PD_global_final.nc', 'r')

PD_global = fh.variables['PD'][:] 	#monthly global area-averaged potential density (kg/m^3)

#Take yearly-averaged rho_ref over first 50 years (or over the SOM cycle?)
PD_global = PD_global[0:50*12,:]
#PD_global = PD_global[(year_start - 1)*12:year_end*12,:]

# Precompute month weights
month_days = np.array([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
month_weights = month_days / np.sum(month_days)
month_weights_tiled = np.tile(month_weights, 50)  # Repeat weights for 50 years
#month_weights_tiled = np.tile(month_weights, year_end - year_start + 1) 

# Compute weighted mean for the years
rho_ref = np.sum(PD_global * month_weights_tiled[:, np.newaxis], axis=0) / np.sum(month_weights_tiled)

#-----------------------------------------------------------------------------------------

# Compute vertical gradient of rho_ref
n0 = np.zeros(len(depth))
for depth_i in range(len(depth)):
	#print(depth_i)
	#print(layer[depth_i])
	#print(dens[depth_i])
	#print(dens[depth_i - 1])
	#print(dens[depth_i + 1])
	#if ma.is_masked(dens[depth_i]):
	#	continue
	if depth_i == 0:
		n0[depth_i] = (rho_ref[depth_i] - rho_ref[depth_i+1])/layer[0]
	elif depth_i == len(depth) - 1:
		n0[depth_i] = (rho_ref[depth_i - 1] - rho_ref[depth_i])/layer[depth_i]
	else:
		n0[depth_i] = (rho_ref[depth_i-1] - rho_ref[depth_i+1])/(0.5*layer[depth_i-1] + layer[depth_i] + 0.5*layer[depth_i+1])

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
fh.variables['n0'].long_name 		= 'Vertical density gradient'

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

#-----------------------------------------------------------------------------------------

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
		PD_month = PD_month - rho_ref[:, np.newaxis, np.newaxis] #rho_star
		PDPD_all[pdw_index, :,:,:] = (PD_month**2) * month_weights[month_i] #rho_star^2
        	
	#Start calculating energetics when year_i has at least 5 years
	if year_i >= window - 1:
	
		#Compute mean values in window
		PD2_mean	= np.sum(PDPD_all, axis = 0) / window #divide by windowlength because we have every months 5 times
		
		APE_year = -0.5 * g * (1 / n0[:, np.newaxis, np.newaxis]) * PD2_mean
		
		#Volume integrated
		APE_int[year_i - window + 1] 	= np.sum(APE_year*volume)
		
	#-----------------------------------------------------------------------------------------
	#-----------------------------------------------------------------------------------------

	print('Data is written to file')
	fh = netcdf.Dataset(directory+'Ocean/APE_year_'+str(year_start)+'-'+str(year_end)+'_window_'+str(window)+'_volume_integrated_SO30_TEST_rhoref_0-50_TEST_goedegsw.nc', 'w')

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
		
