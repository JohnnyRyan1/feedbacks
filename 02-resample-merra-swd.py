#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DESCRIPTION

1. Reproject MERRA-2 data to 1 km

"""

#%%

# Import modules
import xarray as xr
from pyresample import kd_tree, geometry
import numpy as np
#from scipy import ndimage
import netCDF4

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
dest = path + 'merra-swd-resample/'

#%%

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6.nc')

# =============================================================================
# # Covert ISMIP to 50 km grid cells
# lons_50km = ndimage.zoom(ismip_1km['lon'].values, 0.02)
# lats_50km = ndimage.zoom(ismip_1km['lat'].values, 0.02)
# x_50km = ndimage.zoom(ismip_1km['x'].values, 0.02)
# y_50km = ndimage.zoom(ismip_1km['y'].values, 0.02)
# =============================================================================

# Define target projection
target_def = geometry.SwathDefinition(lons=ismip_1km['lon'].values, lats=ismip_1km['lat'].values)

# Define years
years = np.arange(2002, 2024)

#%%

"""
Resample MERRA-2

"""

for year in years:
    
    # Read files
    merra_swd = xr.open_dataset(path + 'merra-swd/swd_'  + str(year) + '.nc')
    
    # Extract lat/lon
    merra_lon, merra_lat = np.meshgrid(merra_swd['longitude'].values, merra_swd['latitude'].values)
    
    # Define source projection
    source_def = geometry.SwathDefinition(lons=merra_lon, lats=merra_lat)
        
    # Reproject
    swd_allsky = kd_tree.resample_nearest(source_def, 
                                   np.rollaxis(merra_swd['swd_allsky'].values, 0, 3), 
                                   target_def, 
                                   radius_of_influence=50000)
    
    # Convert zeros to nans
    swd_allsky[swd_allsky == 0] = np.nan
    
    # Reproject
    swd_clrsky = kd_tree.resample_nearest(source_def, 
                                   np.rollaxis(merra_swd['swd_clrsky'].values, 0, 3), 
                                   target_def, 
                                   radius_of_influence=50000)
    
    # Convert zeros to NaNs
    swd_clrsky[swd_clrsky == 0] = np.nan

    
    ###############################################################################
    # Save 1 km dataset to NetCDF
    ###############################################################################
    dataset = netCDF4.Dataset(dest + 'swd_' + str(year) + '.nc', 
                              'w', format='NETCDF4_CLASSIC')
    print('Creating %s' %dest + 'swd_' + str(year) + '.nc')
    dataset.Title = "Daily downward allsky and clearsky shortwave radiation for %s from MERRA-2" %(str(year))
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C. et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', swd_clrsky.shape[0])
    lon_dim = dataset.createDimension('x', swd_clrsky.shape[1])
    data_dim = dataset.createDimension('z', swd_clrsky.shape[2])

        
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
    allsky_swd_nc = dataset.createVariable('swd_allsky', np.int8, ('y','x','z'))
    clrsky_swd_nc = dataset.createVariable('swd_clrsky', np.int8, ('y','x','z'))

    # Write data to layers
    Y[:] = ismip_1km['lat'].values
    X[:] = ismip_1km['lon'].values
    y[:] = ismip_1km['lat'].values[:,0]
    x[:] = ismip_1km['lon'].values[0,:]
    allsky_swd_nc[:] = swd_allsky.astype(np.float32)
    clrsky_swd_nc[:] = swd_clrsky.astype(np.float32)
    z[:] = np.arange(1,93)
    
    print('Writing data to %s' %dest + 'swd_' + str(year) + '.nc')
        
    # Close dataset
    dataset.close()








