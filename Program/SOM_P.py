#Program determines the SOM-P index

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
	
	#print(lon)
	
	#plt.figure()
	#plt.plot(lon, lon)
	#plt.show()
	
	#Negative longitudes are also allowed
	lon_min 	= 170
	lon_max	        = 250
	lat_min		= -60
	lat_max		= -45
	
	#First get the latitude extent
	lat_min_index	= (fabs(lat - lat_min)).argmin()
	lat_max_index	= (fabs(lat - lat_max)).argmin() + 1	
	lat		= lat[lat_min_index:lat_max_index]
	depth_t		= depth_t[lat_min_index:lat_max_index]
	area		= area[lat_min_index:lat_max_index]
		
	fh = netcdf.Dataset(filename, 'r')

	temp 		= fh.variables['TEMP'][0, lat_min_index:lat_max_index]	#Sea surface temperature

	fh.close()

	#Mask the land
	temp		= ma.masked_where(depth_t <= 0.0, temp)
	
	#Convert the field to 0E to 360E
	lon_2, area	=  ConverterField2D(lon, area)
	lon, temp	= ConverterField2D(lon, temp)
	
	print(lon)

	if lon_min < 0 or lon_max > 360:
		#If needed, extend the field to add negative longitudes and beyond 360E
		lon_2, temp	= PeriodicBoundaries2D(lon, lat, temp, 1800)
		lon, area	= PeriodicBoundaries2D(lon, lat, area, 1800)	
	
	print(lon)
	
	#Get the longitude extent
	lon_min_index	= (fabs(lon - lon_min)).argmin()
	lon_max_index	= (fabs(lon - lon_max)).argmin() + 1	
	lon		= lon[lon_min_index:lon_max_index]
	area		= area[:, lon_min_index:lon_max_index]
	temp		= temp[:, lon_min_index:lon_max_index]

	#Mask the land for the area grid
	area		= ma.masked_array(area, mask = temp.mask)
	
	print(lon)
	print(lat)
	
	plt.figure()
	plt.contourf(lon, lat, temp)
	plt.show()
	
	#sys.exit()

	return lon, lat, area, temp
	
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

lon_min 	= 175
lon_max	        = 210
lat_min		= -60
lat_max		= -45

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/SOM_P_lat_'+str(lat_min)+'-'+str(lat_max)+'_lon_'+str(lon_min)+'-'+str(lon_max)+'.nc', 'w')

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

print('Data is written to file')
