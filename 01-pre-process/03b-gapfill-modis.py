#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Interpolate missing data gaps in MODIS albedo data.

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import glob
import os

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Import MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo-filled/*.nc'))

# Define destination
dest = path + 'modis-albedo-complete/'

# Define mask
mask = ismip_1km['GIMP'].values

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%

pre = []
post = []

for file in modis_files:

    # Get path and filename seperately 
    infilepath, infilename = os.path.split(file) 
    # Get file name without extension            
    infileshortname, extension = os.path.splitext(infilename)
    print('Processing... %s' % infileshortname)
# =============================================================================
#     if os.path.exists(dest + infileshortname + '.nc'):
#         print('Skipping... %s' %(infileshortname + '.nc'))
#     else:
#         print('Processing... %s' % infileshortname)
#         
# =============================================================================
    # Read MODIS data
    f = xr.open_dataset(file)
    lats = f['latitude'].values
    lons = f['longitude'].values
    
    # Convert to float
    array = f['filled_albedo'].values.astype(np.float32)
    
    # Replace zeros with NaNs
    array[array == 0] = np.nan
    
    # Find pixels that contain more than one NaN value
    null_pixels = np.sum(np.isnan(array), axis=2)
    null_pixels = ((mask == 1) & (null_pixels > 0))
    pre.append(np.sum((np.isnan(array)) & (mask3d == 1)))
    
    for j in range(np.where(null_pixels)[0].shape[0]):
        # Identify pixel
        coords = [np.where(null_pixels)[0][j], np.where(null_pixels)[1][j]]
        
        # Define pixel
        pixel = array[coords[0], coords[1], :]
        
        # Find neighbors
        neighbors = array[coords[0]-1:coords[0]+2,coords[1]-1:coords[1]+2, :]
        
        new_values = []
        for i in range(pixel.shape[0]):
            if np.isnan(pixel[i]):
                new_values.append(np.nanmedian(neighbors[:,:,i]))
            else:
                new_values.append(pixel[i])
        
        # Replace values
        array[coords[0], coords[1], :] = new_values
    
    
    # Find pixels that contain more than one NaN value
    null_pixels = np.sum(np.isnan(array), axis=2)
    null_pixels = ((mask == 1) & (null_pixels > 0))
    post.append(np.sum((np.isnan(array)) & (mask3d == 1)))
        
# =============================================================================
#         ###############################################################################
#         # Save 1 km dataset to NetCDF
#         ###############################################################################
#            
#         dataset = netCDF4.Dataset(dest + infileshortname + '.nc', 'w', format='NETCDF4_CLASSIC')
#         print('Creating %s' % dest + infileshortname + '.nc')
#         dataset.Title = "Daily gap-filled albedo for %s from MOD10A1 product" %(str(infileshortname[-4:]))
#         import time
#         dataset.History = "Created " + time.ctime(time.time())
#         dataset.Projection = "WGS 84"
#         dataset.Reference = "Ryan, J. C. et al. (unpublished)"
#         dataset.Contact = "jryan4@uoregon.edu"
#             
#         # Create new dimensions
#         lat_dim = dataset.createDimension('y', array.shape[0])
#         lon_dim = dataset.createDimension('x', array.shape[1])
#         data_dim = dataset.createDimension('z', array.shape[2])
#         
#             
#         # Define variable types
#         Y = dataset.createVariable('latitude', np.float32, ('y','x'))
#         X = dataset.createVariable('longitude', np.float32, ('y','x'))
#         
#         y = dataset.createVariable('y', np.float32, ('y'))
#         x = dataset.createVariable('x', np.float32, ('x'))
#         z = dataset.createVariable('z', np.float32, ('z'))
#             
#         # Define units
#         Y.units = "degrees"
#         X.units = "degrees"
#            
#         # Create the actual 3D variable
#         albedo_nc = dataset.createVariable('albedo', np.int8, ('y','x','z'))
#         
#         # Write data to layers
#         Y[:] = lats
#         X[:] = lons
#         x[:] = lons[0,:]
#         y[:] = lats[:,0]
#         albedo_nc[:] = array.astype(np.int8)
#         z[:] = np.arange(1,93)
#         
#         print('Writing data to %s' % dest + infileshortname + '.nc')
#             
#         # Close dataset
#         dataset.close()
# =============================================================================
    
#%%

# Save DataFrame
df = pd.DataFrame((pre, post)).T
df.columns = ['pre', 'post']

df.to_csv(path + 'filling-stats/gap-filled-totals.csv')





























    
    
    
    
    
    
