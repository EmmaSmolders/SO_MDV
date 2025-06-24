#Program determines the SOM index

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

def ReadinData(filename):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#Grid is rectangular up to 28N
	lon 		= fh.variables['TLONG'][400]		#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['TLAT'][:, 780]		#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	depth_t		= fh.variables['HT'][:] / 100		#Depth of bathymetry
	area		= fh.variables['TAREA'][:] / 10000	#TAREA (m^2)

	
	fh.close()
	
	#Use negative longitudes for 180W - 0W
	lon[lon>180]	= lon[lon>180] - 360.0

	#SOM region, following Le Bars et al. (2016)
	lon_min_index	= (fabs(lon - -50)).argmin()
	lon_max_index	= (fabs(lon - 0)).argmin() + 1
	lat_min_index	= (fabs(lat - -50)).argmin()
	lat_max_index	= (fabs(lat - -35)).argmin() + 1	
	
	lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_t		= depth_t[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	
	fh = netcdf.Dataset(filename, 'r')

	temp 		= fh.variables['TEMP'][0, lat_min_index:lat_max_index, lon_min_index:lon_max_index]	#Sea surface temperature (shape [depth, lat, lon])

	fh.close()

	#Mask the land and the also the area
	temp		= ma.masked_where(depth_t <= 0.0, temp)
	area		= ma.masked_array(area, mask = temp.mask)

	return lon, lat, area, temp
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

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
lon, lat, area, temp	= ReadinData(files[0])
time_year		= ma.masked_all(year_end-year_start+1)
SOM			= ma.masked_all(len(time_year))

#Normalise the area
area_norm		= area / np.sum(area)

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	temp_year 		= ma.masked_all((12, len(lat), len(lon)))
	
	for month_i in range(12):

		#Get the monthly files 
		filename 	= files[year_i*12 + month_i]

#		print(filename)
		lon, lat, area, temp_year[month_i] = ReadinData(filename)

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(lat), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	temp_year	= np.sum(temp_year * month_days_all, axis = 0)
	SOM[year_i]	= np.sum(temp_year * area_norm)
	
	
plot(time_year, SOM)
show()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/SOM.nc', 'w')

fh.createDimension('time', len(time_year))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('SOM', float, ('time'), zlib=True)

fh.variables['SOM'].longname 		= 'The Southern Ocean Mode (SOM)'

fh.variables['time'].units 		= 'year'
fh.variables['SOM'].units 		= 'deg C'

#Writing data to correct variable
fh.variables['time'][:] 		= time_year
fh.variables['SOM'][:] 			= SOM

fh.close()
