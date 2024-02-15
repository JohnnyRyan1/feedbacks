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
import netCDF4

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
dest = path + 'merra-t2m-resample/'

#%%

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6.nc')

# Define target projection
target_def = geometry.SwathDefinition(lons=ismip_1km['lon'].values, lats=ismip_1km['lat'].values)

# Define years
years = np.arange(2002, 2024)

#%%

"""
Resample MERRA-2 SWD

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
    allsky_swd_nc = dataset.createVariable('swd_allsky', np.float32, ('y','x','z'))
    clrsky_swd_nc = dataset.createVariable('swd_clrsky', np.float32, ('y','x','z'))

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

#%%


"""
Resample MERRA-2 surface air temperature

"""

for year in years:
    
    # Read files
    merra_t2m = xr.open_dataset(path + 'merra-t2m/t2m_'  + str(year) + '.nc')
    
    # Extract lat/lon
    merra_lon, merra_lat = np.meshgrid(merra_t2m['longitude'].values, merra_t2m['latitude'].values)
    
    # Define source projection
    source_def = geometry.SwathDefinition(lons=merra_lon, lats=merra_lat)
        
    # Reproject
    t2m = kd_tree.resample_nearest(source_def, 
                                   np.rollaxis(merra_t2m['t2m'].values, 0, 3), 
                                   target_def, 
                                   radius_of_influence=50000)
    
    # Convert zeros to nans
    t2m[t2m == 0] = np.nan

    
    ###############################################################################
    # Save 1 km dataset to NetCDF
    ###############################################################################
    dataset = netCDF4.Dataset(dest + 't2m_' + str(year) + '.nc', 
                              'w', format='NETCDF4_CLASSIC')
    print('Creating %s' %dest + 't2m' + str(year) + '.nc')
    dataset.Title = "Daily air temperature at 500 hPa for %s from MERRA-2" %(str(year))
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C. et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', t2m.shape[0])
    lon_dim = dataset.createDimension('x', t2m.shape[1])
    data_dim = dataset.createDimension('z', t2m.shape[2])

        
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
    t2m_nc = dataset.createVariable('t2m', np.float32, ('y','x','z'))

    # Write data to layers
    Y[:] = ismip_1km['lat'].values
    X[:] = ismip_1km['lon'].values
    y[:] = ismip_1km['lat'].values[:,0]
    x[:] = ismip_1km['lon'].values[0,:]
    t2m_nc[:] = t2m.astype(np.float32)
    z[:] = np.arange(1,93)
    
    print('Writing data to %s' %dest + 't2m_' + str(year) + '.nc')
        
    # Close dataset
    dataset.close()







