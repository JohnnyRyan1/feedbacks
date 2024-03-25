#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Correct clear-sky albedo to all-sky using Key et al. (2001) approach/

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import glob
import os


#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
alt_path = '/Users/' + user + '/Dropbox (University of Oregon)/published/clouds/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo-complete/*.nc'))

# Define MODIS files
other_modis_files = sorted(glob.glob(path + 'modis-albedo/*.nc'))

# Define years
years = np.arange(2002, 2024, 1)

# Define destination to save
dest = path + 'modis-albedo-final/'

#%%

# Produce cloud fraction product
cloud_fraction_data = np.zeros((ismip_1km['GIMP'].shape))

for f in other_modis_files:
    print(f)
    # Import MODIS data
    mod = xr.open_dataset(f)
    
    # Compute cloud fraction
    cloud_frac = (mod['albedo'][:,:,18:110].values == 0).sum(axis=2) / 92
    
    # Add to data
    cloud_fraction_data = np.dstack((cloud_fraction_data, cloud_frac))
    
# Average
cloud_fraction = np.nanmean(cloud_fraction_data[:,:,1:], axis=2)

# Remove values outside of mask
cloud_fraction[ismip_1km['GIMP'] == 0] = np.nan
    
#%%

def save2netcdf(dest, file, lats, lons, albedo, year):
    
    ###############################################################################
    # Save 1 km dataset to NetCDF
    ###############################################################################
    
    # Get path and filename seperately 
    infilepath, infilename = os.path.split(file) 
    # Get file name without extension            
    infileshortname, extension = os.path.splitext(infilename)
        
    dataset = netCDF4.Dataset(dest + infilename, 'w', format='NETCDF4_CLASSIC')
    print('Creating %s' % dest + infilename)
    dataset.Title = "Daily gap-filled, clear-sky corrected albedo for %s from MOD10A1 product" %(str(year))
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C. et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', albedo.shape[0])
    lon_dim = dataset.createDimension('x', albedo.shape[1])
    data_dim = dataset.createDimension('z', albedo.shape[2])

        
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
    filled_albedo_nc = dataset.createVariable('albedo', np.int8, ('y','x','z'))

    # Write data to layers
    Y[:] = lats
    X[:] = lons
    x[:] = lons[0,:]
    y[:] = lats[:,0]
    filled_albedo_nc[:] = albedo.astype(np.int8)
    z[:] = np.arange(1,93)
    
    print('Writing data to %s' % dest + infilename)
        
    # Close dataset
    dataset.close()

#%%

for f in modis_files:
    
    # Get path and filename seperately 
    infilepath, infilename = os.path.split(f) 
    # Get file name without extension            
    infileshortname, extension = os.path.splitext(infilename)
    print('Processing... %s' % f)

    # Import MODIS data
    mod = xr.open_dataset(f)
    
    # Correct
    corrected = mod['albedo'].values + (cloud_fraction[:,:,np.newaxis] * 5)

    # Save as NetCDF
    save2netcdf(dest, f, mod['latitude'].values, mod['longitude'].values, 
                corrected, infileshortname[-4:])





