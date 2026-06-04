#Maximum mixed layer depth of Drake Passage region

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

	fh = netcdf.Dataset("/projects/0/prace_imau/prace_2013081679/pop/tx0.1v2/pop.B2000.tx0.1v2.qe_hosing.001/rcp8.5_co2_f05_t12.pop.h.2001-01.nc", 'r')

	#Grid is rectangular up to 28N (TAUX and TAUY are on U-grid)
	lon 		= fh.variables['TLONG'][400]		#Longitude (at lat-index 400 there is no land mask, continious longitudes)
	lat 		= fh.variables['TLAT'][:, 780]		#Latitude  (at lon-index 780 there is no land mask, continious latitudes)
	depth_t		= fh.variables['HT'][:] / 100		#Depth of bathymetry
	area		= fh.variables['TAREA'][:] / 10000	#TAREA (m^2)

	fh.close()
	
	#plt.figure()
	#plt.plot(lon, lon)
	#plt.show()
	
	#Use negative longitudes for 180W - 0W
	#lon[lon>180]	= lon[lon>180] - 360.0

	#Select region (AU region)
	#lon_min_index	= (fabs(lon - 80)).argmin()
	#lon_max_index	= (fabs(lon - 150)).argmin() + 1
	#lat_min_index	= (fabs(lat - -70)).argmin()
	#lat_max_index	= (fabs(lat - -50)).argmin() + 1	
	
	#lon		= lon[lon_min_index:lon_max_index]
	#lat		= lat[lat_min_index:lat_max_index]
	#depth_t		= depth_t[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	#area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	
	#Select region (NZ region) Swithc off using negative longitudes
	lon_min_index	= (fabs(lon - 160)).argmin()
	lon_max_index	= (fabs(lon - 190)).argmin() + 1
	lat_min_index	= (fabs(lat - -70)).argmin()
	lat_max_index	= (fabs(lat - -60)).argmin() + 1	
	
	lon		= lon[lon_min_index:lon_max_index]
	lat		= lat[lat_min_index:lat_max_index]
	depth_t		= depth_t[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	area		= area[lat_min_index:lat_max_index, lon_min_index:lon_max_index]
	
	#plt.figure()
	#plt.plot(lon, lon)
	#plt.show()
	
	#sys.exit()
	
	fh = netcdf.Dataset(filename, 'r')

	#Get the relevant region
	mld 		= fh.variables['HMXL'][lat_min_index:lat_max_index, lon_min_index:lon_max_index]/100		#Mixed layer depth [m]
#	mld 		= fh.variables['HMXL'][lat_min_index:lat_max_index, :]/100
	fh.close()
	
	#Remove land mask
	mld		= ma.masked_where(depth_t <= 0, mld)

	print(np.shape(lat))
	print(np.shape(lon))
	print(np.shape(mld))

	#plt.figure()
	#plt.contourf(lon, lat, mld)
	#plt.savefig(directory + 'test_mld.pdf')
	#plt.show()

	#sys.exit()

	return lat, lon, area, mld
   			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

#Second SOM cycle (model year 45-91), last SOM cycle (495-571)
year_start	= 1
year_end	= 600

files = []

files_all	= glob.glob(directory_data+'t.t0.1_42l_nccs01.*.nc')
files_all.sort()

print(files_all[0])
print(files_all[-1])
print(len(files_all[0]))

for file_i in range(len(files_all)):

	if len(files_all[file_i]) == 116:
		files.append(files_all[file_i])
		
print(files[0])
print(files[-1]) 
print(year_start*12 - 1)
print(year_end*12)

#note that january 300 does not exist and that february 300 does not have PD data. Therefore, we take january 301 to january 351
files	= files[:600*12]

print(files[0])
print(files[-1])

#-----------------------------------------------------------------------------------------

lat, lon, area, mld		= ReadinData(files[0])
time_year			= ma.masked_all(year_end-year_start+1)
mld_all				= ma.masked_all((len(time_year), len(lat), len(lon)))
mld_max				= ma.masked_all(len(time_year))

print(lat)
print(lon)
print(np.shape(lat))
print(np.shape(lon))
print(np.shape(mld))

#Store maximum MLD in WKGP region
for year_i in range(len(time_year)):
	#Now determine for each month
	print(year_i)
	time_year[year_i] 	= year_i + year_start

	mld_year 		= ma.masked_all((12, len(lat), len(lon)))
	
	for month_i in range(12):

		#Get the monthly files
		filename 	= files[year_i*12 + month_i]

		print(filename)
		lat, lon, area, mld_year[month_i] = ReadinData(filename)

	#------------------------------------------------------------------------------
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(lat), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#-----------------------------------------------------------------------------------------

	#Determine the time mean over the months of choice
	mld_all[year_i]	= np.sum(mld_year * month_days_all, axis = 0)
	
	#Determine maximum MLD in WKGP region
	mld_max[year_i] = np.max(mld_all[year_i,:,:])

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/MLD_year_'+str(year_start)+'-'+str(year_end)+'_NZ_160_190E_70S_62S.nc', 'w')

fh.createDimension('time', len(time_year))
fh.createDimension('lat', len(lat))
fh.createDimension('lon', len(lon))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('MLD_max', float, ('time'), zlib=True)
fh.createVariable('MLD', float, ('time', 'lat', 'lon'), zlib=True)
fh.createVariable('lat', float, ('lat'), zlib=True)
fh.createVariable('lon', float, ('lon'), zlib=True)

fh.variables['time'].long_name 		= 'Model year'
fh.variables['MLD_max'].long_name 	= 'Maximum Mixed Layer Depth in region'
fh.variables['MLD'].long_name		= 'Yearly averaged MLD'
fh.variables['lat'].long_name		= 'Latitudes'
fh.variables['lon'].long_name		= 'Longitudes'

fh.variables['time'].units 		= 'year'
fh.variables['MLD_max'].units 		= 'm'
fh.variables['MLD'].units		= 'm'

#Writing data to correct variable
fh.variables['time'][:] 		= time_year
fh.variables['MLD_max'][:] 		= mld_max
fh.variables['MLD'][:]			= mld_all
fh.variables['lat'][:]			= lat
fh.variables['lon'][:]			= lon

fh.close()













