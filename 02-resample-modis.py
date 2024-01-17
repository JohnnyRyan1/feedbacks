#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Resample MODIS data to 1 km ISMIP grid.

"""

# Import modules
import netCDF4
import numpy as np
import glob
from pyresample import geometry, utils, image
import os

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Define destination to save
dest_1km = path + 'modis-albedo/'

# Import ISMIP 1 km grid
ismip_1km = netCDF4.Dataset(path + '1km-ISMIP6.nc','r')

# Define MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo-intermediate/*.nc'))

# Define years
years = np.arange(2002, 2024, 1)

#%%

for i in years:
    
    if os.path.exists(dest_1km + 'mod10a1-albedo-' + str(i) + '.nc'):
        print('Skipping... %s' %(dest_1km + 'mod10a1-albedo-' + str(i) + '.nc'))
    else:
        modis_list = []   
        # Get MODIS tiles
        for f in modis_files:
            if  f[-7:-3]  == str(i):
                modis_list.append(f)
        
        # Define new master grid
        master_grid_albedo = np.zeros((7200,7200,125), dtype='float')
        master_grid_lat = np.zeros((7200, 7200), dtype='float')
        master_grid_lon = np.zeros((7200, 7200), dtype='float')
    
        # Add tile to master grid
        for j in modis_list:
            if j[-13:-8] == '15v00':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[0:2400,0:2400,:] = modis.variables['albedo'][:]
                master_grid_lat[0:2400,0:2400] = modis.variables['latitude'][:]
                master_grid_lon[0:2400,0:2400] = modis.variables['longitude'][:]
            if j[-13:-8] == '16v00':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[0:2400,2400:4800,:] = modis.variables['albedo'][:]
                master_grid_lat[0:2400,2400:4800] = modis.variables['latitude'][:]
                master_grid_lon[0:2400,2400:4800] = modis.variables['longitude'][:]
            if j[-13:-8] == '17v00':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[0:2400,4800:7200,:] = modis.variables['albedo'][:]
                master_grid_lat[0:2400,4800:7200] = modis.variables['latitude'][:]
                master_grid_lon[0:2400,4800:7200] = modis.variables['longitude'][:]
            if j[-13:-8] == '15v01':
                modis = netCDF4.Dataset(j, 'r')  
                master_grid_albedo[2400:4800,0:2400,:] = modis.variables['albedo'][:]
                master_grid_lat[2400:4800,0:2400] = modis.variables['latitude'][:]
                master_grid_lon[2400:4800,0:2400] = modis.variables['longitude'][:]
            if j[-13:-8] == '16v01':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[2400:4800,2400:4800,:] = modis.variables['albedo'][:]
                master_grid_lat[2400:4800,2400:4800] = modis.variables['latitude'][:]
                master_grid_lon[2400:4800,2400:4800] = modis.variables['longitude'][:]
            if j[-13:-8] == '17v01':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[2400:4800,4800:7200,:] = modis.variables['albedo'][:]
                master_grid_lat[2400:4800,4800:7200] = modis.variables['latitude'][:]
                master_grid_lon[2400:4800,4800:7200] = modis.variables['longitude'][:]
            if j[-13:-8] == '15v02':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[4800:7200,0:2400,:] = modis.variables['albedo'][:]
                master_grid_lat[4800:7200,0:2400] = modis.variables['latitude'][:]
                master_grid_lon[4800:7200,0:2400] = modis.variables['longitude'][:]
            if j[-13:-8] == '16v02':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[4800:7200,2400:4800,:] = modis.variables['albedo'][:]
                master_grid_lat[4800:7200,2400:4800] = modis.variables['latitude'][:]
                master_grid_lon[4800:7200,2400:4800] = modis.variables['longitude'][:]
            if j[-13:-8] == '17v02':
                modis = netCDF4.Dataset(j, 'r')
                master_grid_albedo[4800:7200,4800:7200,:] = modis.variables['albedo'][:]
                master_grid_lat[4800:7200,4800:7200] = modis.variables['latitude'][:]
                master_grid_lon[4800:7200,4800:7200] = modis.variables['longitude'][:]
    
        # Get ISMIP6 lat lons
        lon_1km = ismip_1km.variables['lon'][:]
        lat_1km = ismip_1km.variables['lat'][:]
        
        # Convert 0s to NaNs so they do not interfere with resampling
        master_grid_albedo = master_grid_albedo.astype('float')
        master_grid_albedo[master_grid_albedo == 0] = np.nan
        
        # Define regridding conversion
        swath_def = geometry.SwathDefinition(lons=lon_1km, lats=lat_1km)
        swath_con = geometry.SwathDefinition(lons=master_grid_lon, lats=master_grid_lat)
        albedo_con = image.ImageContainer(master_grid_albedo, swath_con)
        row_indices, col_indices = utils.generate_nearest_neighbour_linesample_arrays(swath_con, swath_def, 1000)
        
        # Perform regridding
        albedo_result = albedo_con.get_array_from_linesample(row_indices, col_indices)
        
        # Filter bad values
        albedo_result[albedo_result < 0] = 0
        albedo_result[albedo_result > 100] = 0
           
        ###############################################################################
        # Save 1 km dataset to NetCDF
        ###############################################################################
        dataset = netCDF4.Dataset(dest_1km + 'mod10a1-albedo-' + str(i) + '.nc', 
                                  'w', format='NETCDF4_CLASSIC')
        print('Creating %s' %dest_1km + 'mod10a1-albedo-' + str(i) + '.nc')
        dataset.Title = "Daily albedo for %s from MOD10A1 product" %(str(i))
        import time
        dataset.History = "Created " + time.ctime(time.time())
        dataset.Projection = "WGS 84"
        dataset.Reference = "Ryan, J. C. et al. (unpublished)"
        dataset.Contact = "jryan4@uoregon.edu"
            
        # Create new dimensions
        lat_dim = dataset.createDimension('y', albedo_result.shape[0])
        lon_dim = dataset.createDimension('x', albedo_result.shape[1])
        data_dim = dataset.createDimension('z', albedo_result.shape[2])
    
            
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
        albedo_nc = dataset.createVariable('albedo', np.int8, ('y','x','z'))
        
        # Write data to layers
        Y[:] = lat_1km
        X[:] = lon_1km
        x[:] = lon_1km[0,:]
        y[:] = lat_1km[:,0]
        albedo_nc[:] = albedo_result.astype(np.int8)
        z[:] = np.arange(1,126)
        
        print('Writing data to %s' %dest_1km + 'mod10a1-albedo-' + str(i) + '.nc')
            
        # Close dataset
        dataset.close()
