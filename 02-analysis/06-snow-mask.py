#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute a snow mask.

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import glob
import itertools
import os
import netCDF4

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo-complete/*.nc'))

#%%

# Define mask
mask = ismip_1km['GIMP'].values

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)


#%%

# Remove exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded_grid_cells.csv')

# Filter
mask[(nan_df['x'].values, nan_df['y'].values)] = False

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%
snow_mask = np.zeros(mask.shape)

for f in range(len(modis_files)):

    print('Processing... %s' % modis_files[f])

    # Import albedo data
    modis = xr.open_dataset(modis_files[f])

    # Some preprocessing
    albedo = modis['albedo'].values.astype(np.float32)
    albedo[mask3d == 0] = np.nan
    albedo[albedo == 0] = np.nan
            
    # Add min and max snow albedo
    albedo[albedo > 83] = 83
    albedo[albedo <= 30] = 30

    # Classify
    classified = np.copy(albedo)
    classified[classified <= 55] = 1
    classified[classified > 55] = 2
    
    # Identify pixels which are always snow
    mean_classified = np.nanmean(classified, axis=2)

    # Stack
    snow_mask = np.dstack((snow_mask, (mean_classified == 2)))

# Identify pixels tha are always snow through the entire study period
snow_mask = snow_mask[:,:,1:]
snow_mask_mean = np.nanmean(snow_mask, axis=2)
snow_mask_snow = (snow_mask_mean == 1)

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
   
lats = modis['latitude'].values
lons = modis['longitude'].values

dataset = netCDF4.Dataset(path + 'snow-mask.nc', 'w', format='NETCDF4_CLASSIC')
print('Creating %s' % path + 'snow-mask.nc')
dataset.Title = "Map of pixels that are always snow during the study period"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', snow_mask_snow.shape[0])
lon_dim = dataset.createDimension('x', snow_mask_snow.shape[1])

    
# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))

y = dataset.createVariable('y', np.float32, ('y'))
x = dataset.createVariable('x', np.float32, ('x'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
snow_mask_nc = dataset.createVariable('snow_mask', np.int8, ('y','x'))

# Write data to layers
Y[:] = lats
X[:] = lons
x[:] = lons[0,:]
y[:] = lats[:,0]
snow_mask_nc[:] = snow_mask_snow.astype(np.int8)

print('Writing data to %s' % path + 'snow-mask.nc')
    
# Close dataset
dataset.close()














