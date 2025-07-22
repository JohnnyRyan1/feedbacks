#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute surface albedo climatology.

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import glob

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo-final/*.nc'))

# Define years
years = np.arange(2002, 2024, 1)

#%%

albedo_climatology = np.zeros((2881, 1681, 1))

# Produce SWD climatology
for i in range(92):
    print(i)
    clim = np.zeros((2881, 1681, 1))
    
    for year in years:
        
        # Read files
        modis = xr.open_dataset(path + 'modis-albedo-final/mod10a1-albedo-'  + str(year) + '.nc')
    
        # Stack 
        clim = np.dstack((clim, modis['albedo'][:,:,i].values))
    
    # Remove first layer
    clim = clim[:,:,1:]
    
    # Stack to main climatology
    albedo_climatology = np.dstack((albedo_climatology, np.nanmean(clim, axis=2)))

# Remove first layer
albedo_climatology = albedo_climatology[:,:,1:]

#%%

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'modis-climatology.nc', 
                          'w', format='NETCDF4_CLASSIC')
print('Creating %s' %path + 'modis-climatology.nc')
dataset.Title = "Climatology of daily albedo for 2002-2023 period from MODIS"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', albedo_climatology.shape[0])
lon_dim = dataset.createDimension('x', albedo_climatology.shape[1])
data_dim = dataset.createDimension('z', albedo_climatology.shape[2])

    
# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))

y = dataset.createVariable('y', np.float32, ('y'))
x = dataset.createVariable('x', np.float32, ('x'))
z = dataset.createVariable('z', np.float32, ('z'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
albedo_nc = dataset.createVariable('albedo_climatology', np.int8, ('y','x','z'))

# Write data to layers
Y[:] = ismip_1km['lat'].values
X[:] = ismip_1km['lon'].values
y[:] = ismip_1km['lat'].values[:,0]
x[:] = ismip_1km['lon'].values[0,:]
albedo_nc[:] = albedo_climatology.astype(np.int8)
z[:] = np.arange(1,93)

print('Writing data to %s' %path + 'modis-climatology.nc')
    
# Close dataset
dataset.close()


#%%













