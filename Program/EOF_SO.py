#Program conduct the EOF analysis for the Southern Ocean sea surface temperatures

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory_data	= '/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/tavg/'
directory	= '/home/smolders/HR_POP/Data/Ocean/'

def ReadinData(filename, lon_min, lon_max, lat_min, lat_max, lon_lat_red = 1):
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

	if lon_min < 0 or lon_max > 360:
		#If needed, extend the field to add negative longitudes and beyond 360E
		lon_2, temp	= PeriodicBoundaries2D(lon, lat, temp, 1800)
		lon, area	= PeriodicBoundaries2D(lon, lat, area, 1800)	
	
	#Get the longitude extent
	lon_min_index	= (fabs(lon - lon_min)).argmin()
	lon_max_index	= (fabs(lon - lon_max)).argmin() + 1	
	lon		= lon[lon_min_index:lon_max_index]
	area		= area[:, lon_min_index:lon_max_index]
	temp		= temp[:, lon_min_index:lon_max_index]

	#Mask the land for the area grid
	area		= ma.masked_array(area, mask = temp.mask)
	
	if lon_lat_red > 1:
		#Reduce the lon-lat grid by the given factor
		lon_red		= np.zeros(int(len(lon) / lon_lat_red))
		lat_red		= np.zeros(int(len(lat) / lon_lat_red))
		area_red	= ma.masked_all((len(lat_red), len(lon_red)))
		temp_red	= ma.masked_all((len(lat_red), len(lon_red)))

		for lat_i in range(len(lat_red)):
			for lon_i in range(len(lon_red)):
				#Loop over the reduced longitude and latitude grid
				
				#Get the area of the corresponding grids, and normalise
				area_grid	= area[lat_i*lon_lat_red:(lat_i+1)*lon_lat_red, lon_i*lon_lat_red:(lon_i+1)*lon_lat_red]
				area_grid_norm	= area_grid / np.sum(area_grid)
				
				#Take the normalised temperature
				temp_red[lat_i, lon_i]	= np.sum(temp[lat_i*lon_lat_red:(lat_i+1)*lon_lat_red, lon_i*lon_lat_red:(lon_i+1)*lon_lat_red] * area_grid_norm)
				
				#Take the area sum, such that the reduced grid have a consistent area
				area_red[lat_i, lon_i]	= np.sum(area_grid)
				
				#For the coordinates, simply the mean can be taken, as the grid is rectangular
				lon_red[lon_i]		= np.mean(lon[lon_i*lon_lat_red:(lon_i+1)*lon_lat_red])
				lat_red[lat_i]		= np.mean(lat[lat_i*lon_lat_red:(lat_i+1)*lon_lat_red])

		#Return the reduced arrays instead
		return lon_red, lat_red, area_red, temp_red
		
	#plt.figure()
	#plt.contourf(lon, lat, temp)
	#plt.show()
	
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

#lon_min		= 80	#Negative longitudes are allowed
#lon_max		= 150
#lon_min     = -50
#lon_max     = 50	
#lon_min     = 170
#lon_max     = 300
lon_min = 0
lon_max = 1000
lat_min		= -70
lat_max		= -35

lon_lat_red	= 10	#Reduce the longitude and latitude by the given factor, when 1 = no reduction

year_start	= 1
year_end	= 100

month_start 	= 1
month_end 	= 12

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])

files	= files[:600*12]
files	= files[(year_start-1)*12:year_end*12]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------

lon, lat, area, temp	= ReadinData(files[0], lon_min, lon_max, lat_min, lat_max, lon_lat_red)
time_year		= ma.masked_all(year_end-year_start+1)
temp_all		= ma.masked_all((len(time_year), len(lat), len(lon)))	

print(len(lon))
print(lon[0])
print(lon[-1])

#plt.figure()
#plt.contourf(lon, lat, temp)
#plt.show()

#sys.exit()

for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start
	temp_year 		= ma.masked_all((12, len(lat), len(lon)))
	
	for month_i in range(12):
		#Get the monthly files
		filename	= files[year_i*12+month_i]

		print(filename)
		lon, lat, area, temp_year[month_i] = ReadinData(filename, lon_min, lon_max, lat_min, lat_max, lon_lat_red)

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(lat), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	temp_all[year_i]	= np.sum(temp_year * month_days_all, axis = 0)
	
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Please insert your EOF script here <3
#Don't forgot to make your data ready before performing the EOF analysis
#Such as removing trends, normalising the grids to their area etc.

#%% Functions

#Central moving average
def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def TrendRemover(time, data, trend_type):
	"""Removes trend of choice"""
	
	rank = polyfit(time, data, trend_type)
	fitting = 0.0 
		
	for rank_i in range(len(rank)):
			
		fitting += rank[rank_i] * (time**(len(rank) - 1 - rank_i))

	data -= fitting
	
	return data

def MonthRemover(time, data_all):
    """Removes the monthly average for each time series"""
    
    for month_i in range(12):
        month_index             = np.arange(month_i, len(time), 12)
        data_all[month_index]   = data_all[month_index] - np.mean(data_all[month_index], axis = 0)

    return data_all

def MovingAverage(time, data, moving_average):
	"""Determines moving average of time series"""
	
	#Empty array's for the smoothed time series
	time_2     = ma.masked_all(len(time) - moving_average + 1)
	data_2 	   = ma.masked_all(len(time_2))
	
	#Determine the so-called middle index where the moving average is determined
	selection  = (moving_average - 1) / 2		

	for time_i in range(selection, selection + len(time_2)):
		#Take moving average
		data_2[time_i - selection]  = np.mean(data[time_i-selection:time_i-selection + moving_average], axis = 0)
		time_2[time_i - selection] 	= time[time_i]
	
	return time_2, data_2

def Distance(lon1, lat1, lon2, lat2):
	"""Returns distance (m) of two points located at the globe
	coordinates need input in degrees"""

	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2]) #Convert to radians

	#Haversine formula 
	dlon = lon2 - lon1 
	dlat = lat2 - lat1 
	a = math.sin(dlat/2.0)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2.0)**2
	c = 2.0 * math.asin(sqrt(a)) 
	r = 6371000.0 # Radius of earth in meters
	
	return c * r #Distance between two points in meter

moving_average = 5

# Define a function for EOF analysis
def perform_eof_analysis(SST, time, lat, lon, area, trend_type=1, remove_month=0, moving_average=moving_average, eigen_number=1, total_variance=90):
    """
    Perform EOF analysis on the given SST dataset.
    """
    
    #plt.figure()
    #plt.contourf(lon, lat, SST[0,:,:])
    #plt.show()
    
    # Area-weighted average
    area_grids = area
    area_grids = ma.masked_array(area_grids, mask=SST[0].mask)
    area_grids = area_grids / area_grids.max()  # Normalize area

    # Pre-processing
    if trend_type > 0:
        print('Trend is removed')
        for lat_i in range(len(lat)):
            for lon_i in range(len(lon)):
                SST[:, lat_i, lon_i] = TrendRemover(time, SST[:, lat_i, lon_i], trend_type)

    if remove_month == 1:
        print('Monthly signal is removed')
        SST = MonthRemover(time, SST)

    if moving_average > 1:
        print('Moving average is applied')
        time_2 = np.zeros((len(time) - moving_average + 1))
        SST_2 = ma.masked_all((len(time_2), len(lat), len(lon)))

        selection = int((moving_average - 1) / 2)

        for time_i in range(int(selection), int(selection + len(time_2))):
            time_2[time_i - selection] = time[time_i]
            SST_2[time_i - selection] = np.mean(SST[time_i - selection:time_i - selection + moving_average], axis=0)

        time = time_2
        SST = SST_2
        del time_2, SST_2

    print('Data is normalized')
    SST = SST - np.mean(SST, axis=0)
    SST = SST / np.std(SST, axis=0)

    print('Data is scaled by area')
    SST = SST * area_grids

    # EOF analysis
    masked_field = SST[0].mask
    data_all = ma.masked_all((len(lon) * len(lat) - sum(masked_field), len(time)))
    grid_counter = 0

    for lat_i in range(len(lat)):
        for lon_i in range(len(lon)):
            if masked_field[lat_i, lon_i] == False:
                data_all[grid_counter] = SST[:, lat_i, lon_i]
                grid_counter += 1

    print('Determining EOFs and PCs')
    u1, s1, v1 = np.linalg.svd(data_all)
    
    print(s1)
    print(u1)
    print(v1)

    # Determine the variance
    s1 = s1 ** 2.0
    print(eigen_number)
    print(s1)
    print('Eigenvector number', eigen_number, 'contributes ', round(s1[eigen_number - 1] / sum(s1) * 100.0, 1), '% of the total variance')

    variance = 0.0
    for eigen_i in range(len(u1)):
        variance += s1[eigen_i] / sum(s1) * 100.0
        if round(variance, 1) >= total_variance:
            print(total_variance)
            print(variance)
            print(len(u1))
            print('Number of EOFs needed to include at least ' + str(total_variance) + '% of the variance:', eigen_i, '')
            break

    # Retrieve the EOFs
    eof = ma.masked_all((5, len(lat), len(lon)))
    grid_counter = 0

    for lat_i in range(len(lat)):
        for lon_i in range(len(lon)):
            if masked_field[lat_i, lon_i] == False:
                eof[:, lat_i, lon_i] = u1[grid_counter, :len(eof)]
                grid_counter += 1

    return eof, u1, s1, v1, time

#%%
# Process all SST datasets
#lon1, lon2 = -90, 30
#lat1, lat2 = 0, 80
#time1, time2 = 0, 2200
#depth_level = 500

#lat_min_index  = (np.abs(lat[:,0] - lat1)).argmin()
#lat_max_index  = (np.abs(lat[:,0] - lat2)).argmin()+1
#lon_min_index   = (np.abs(lon[0,:] - lon1)).argmin()
#lon_max_index   = (np.abs(lon[0,:] - lon2)).argmin()+1

datasets = {
    "SST": (temp_all, time_year)}

results = {}

for name, (SST, time) in datasets.items():
    print(f"Processing {name}...")
    print(np.shape(SST))
    print(np.shape(time))
    print(np.shape(lat))
    print(np.shape(lon))
    print(np.shape(area))
    eof, u1, s1, v1, time = perform_eof_analysis(SST, time, lat, lon, area)
    results[name] = {
        "eof": eof,
        "u1": u1,
        "s1": s1,
        "v1": v1,
    }
    print(f"Finished processing {name}.")
    print(np.shape(v1))
    print(np.shape(eof))
    print(v1)
    print(np.nanmean(eof))

    # Save the results to a NetCDF file
    filename = f"{directory}EOF_SOM_{name}_month_{month_start}_{month_end}_moving_average_{moving_average}_POP_year_{int(time[0])}_{int(time[-1])}_Pacific.nc"
    print(f"Saving results to {filename}...")

    fh = netcdf.Dataset(filename, 'w')
    fh.createDimension('lon', len(lon))
    fh.createDimension('lat', len(lat))
    fh.createDimension('eof', len(eof))
    fh.createDimension('time', len(time))

    fh.createVariable('lon', float, ('lon'), zlib=True)
    fh.createVariable('lat', float, ('lat'), zlib=True)
    fh.createVariable('eof', float, ('eof'), zlib=True)
    fh.createVariable('time', float, ('time'), zlib=True)
    fh.createVariable('PC', float, ('eof', 'time'), zlib=True)
    fh.createVariable('VAR', float, ('eof'), zlib=True)
    fh.createVariable('EOF', float, ('eof', 'lat', 'lon'), zlib=True)

    fh.variables['lon'].longname = 'Array of longitudes'
    fh.variables['lat'].longname = 'Array of latitudes'
    fh.variables['PC'].long_name = 'PCs of the salinity'
    fh.variables['VAR'].long_name = 'Variance of the PCs/EOFS'
    fh.variables['EOF'].long_name = 'EOFs of the salinity'

    fh.variables['time'].units = 'Model year'
    fh.variables['lon'].units = 'Degrees east'
    fh.variables['lat'].units = 'Degrees north'
    fh.variables['VAR'].units = '%'

    # Writing data to correct variable
    fh.variables['lon'][:] = lon
    fh.variables['lat'][:] = lat
    fh.variables['time'][:] = time
    fh.variables['eof'][:] = np.arange(len(eof)) + 1
    fh.variables['PC'][:] = v1[:len(eof)]
    fh.variables['VAR'][:] = s1[:len(eof)] / sum(s1) * 100.0
    fh.variables['EOF'][:] = eof

    fh.close()
    print(f"Results saved to {filename}.\n")




