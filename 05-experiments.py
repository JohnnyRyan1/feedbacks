#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute SWnet for five experiments.

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import glob
import itertools
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

# Define MERRA files
merra_files = sorted(glob.glob(path + 'merra-swd-resample/*.nc'))

#%%

# Define ice threshold
ice_threshold = [53, 55, 57]

# Define albedo uncertainty
uncertainty = [-2, 0, 2]

#%%

# Define mask
mask = ismip_1km['GIMP'].values

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

# Define ice threshold
ice_threshold = [55]

# Define albedo uncertainty
uncertainty = [0]

#%%

"""
OPTIONAL
"""

# Do some preprocessing to find which grid cells to mask
nan_x, nan_y = [], []
for f in range(len(modis_files)):
    
    print('Processing... %s' % modis_files[f])

    # Import albedo data
    modis = xr.open_dataset(modis_files[f])

    # Some preprocessing
    albedo = modis['albedo'].values.astype(np.float32)
    albedo[mask3d == 0] = np.nan
    albedo[albedo == 0] = np.nan

    # Find NaN values within ice sheet mask
    nan_values = np.where(np.isnan(np.mean(albedo, axis=2)) & (mask == 1))

    # Append
    nan_x.append(nan_values[0])
    nan_y.append(nan_values[1])

# Make into list
flat_nan_x = np.array([num for elem in nan_x for num in elem])
flat_nan_y = np.array([num for elem in nan_y for num in elem])

# Update mask
mask[(flat_nan_x, flat_nan_y)] = False

# Remove duplicates to derive true count
nan_df = pd.DataFrame((flat_nan_x, flat_nan_y)).T
nan_df.columns = ['x', 'y']

# Group
nan_df_group = nan_df.groupby(['x','y']).size().reset_index(name='count')
print('Number of grid cells excluded from analysis = %s' % str(nan_df_group.shape[0]))

# Save as DataFrame
nan_df_group.to_csv(path + 'excluded_grid_cells.csv')

#%%

# Remove exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded_grid_cells.csv')

# Filter
mask[(nan_df['x'].values, nan_df['y'].values)] = False

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%

exp1, exp2, exp3, exp4, exp5 = [], [], [], [], []
ice_diff, snow_diff, snowline_diff = np.zeros(mask.shape), np.zeros(mask.shape), np.zeros(mask.shape)
bulk_diff, bulk_swnet = np.zeros(mask.shape), np.zeros(mask.shape)

for i, j in itertools.product(ice_threshold, uncertainty):

    for f in range(len(modis_files)):

        print('Processing... %s' % modis_files[f])

        # Import albedo data
        modis = xr.open_dataset(modis_files[f])

        # Import SW data
        merra = xr.open_dataset(merra_files[f])

        # Some preprocessing
        albedo = modis['albedo'].values.astype(np.float32)
        albedo[mask3d == 0] = np.nan
        albedo[albedo == 0] = np.nan
                
        # Add max snow albedo
        albedo[albedo > 83] = 83
        albedo[albedo <= 30] = 30

        # Classify
        classified = np.copy(albedo)
        classified[classified <= i] = 1
        classified[classified > i] = 2

        #######################################################################
        # Observed albedo
        #######################################################################
        observed = np.copy(albedo)
        observed = observed + j
    
        # Compute SWnet
        swnet = (1 - (observed / 100)) * merra['swd_allsky']
        swnet_mean = np.nanmean(swnet.values, axis=2)

        # Mask ice sheet
        swnet_mean[mask == 0] = np.nan
        
        #######################################################################
        # Fix snow albedo
        #######################################################################
        
        # Fix snow albedo to fresh snow
        fix_snow = np.copy(albedo)
        #june1 = np.repeat(fix_snow[:,:,0:1], 92, axis=2)
        
        
        fix_snow[classified == 2] = 83 + j
        fix_snow[classified == 1] = fix_snow[classified == 1] + j
        
        # Compute SWnet
        swnet_fixed_snow = (1 - (fix_snow / 100)) * merra['swd_allsky']
        swnet_fixed_snow_mean = np.nanmean(swnet_fixed_snow.values, axis=2)
                                
        # Mask ice sheet
        swnet_fixed_snow_mean[mask == 0] = np.nan
        
        #######################################################################
        # Fix glacier ice albedo
        #######################################################################
        
        # Fix glacier ice albedo to the ice threshold
        fix_ice = np.copy(albedo)
        fix_ice[classified == 1] = i
        fix_ice[classified == 2] = fix_ice[classified == 2] + j
        
        # Compute SWnet
        swnet_fixed_ice = (1 - (fix_ice / 100)) * merra['swd_allsky']
        swnet_fixed_ice_mean = np.nanmean(swnet_fixed_ice.values, axis=2)
                                
        # Mask ice sheet
        swnet_fixed_ice_mean[mask == 0] = np.nan
        
        #######################################################################
        # Fix snow albedo to that 83 and glacier ice albedo to 55
        #######################################################################
        
        # Fix both
        fix_both = np.copy(albedo)
        fix_both[classified == 1] = i
        fix_both[classified == 2] = 83 + j
        
        # Compute SWnet
        swnet_fixed_both = (1 - (fix_both / 100)) * merra['swd_allsky']
        swnet_fixed_both_mean = np.nanmean(swnet_fixed_both.values, axis=2)
                                
        # Mask ice sheet
        swnet_fixed_both_mean[mask == 0] = np.nan
        
        #######################################################################
        # Fix all the ice sheet to 83
        #######################################################################
        
        # Fix all with snow and glacier ice albedo observed on June 1
        fix_all = np.copy(albedo)
        fix_all[classified > 0] = 83 + j
        
        # Compute SWnet
        swnet_fixed_all = (1 - (fix_all / 100)) * merra['swd_allsky']
        swnet_fixed_all_mean = np.nanmean(swnet_fixed_all.values, axis=2)
        
        # Mask ice sheet
        swnet_fixed_all_mean[mask == 0] = np.nan
        
        # Append to lists
        exp1.append(np.nansum(swnet_fixed_all_mean))
        exp2.append(np.nansum(swnet_mean))
        exp3.append(np.nansum(swnet_fixed_both_mean))
        exp4.append(np.nansum(swnet_fixed_ice_mean))
        exp5.append(np.nansum(swnet_fixed_snow_mean))
        
        # Append to grids
        ice_diff = np.dstack((ice_diff, swnet_mean - swnet_fixed_ice_mean))
        snow_diff = np.dstack((snow_diff, swnet_mean - swnet_fixed_snow_mean))
        snowline_diff = np.dstack((snowline_diff, swnet_fixed_both_mean - swnet_fixed_all_mean))
        bulk_diff = np.dstack((bulk_diff, swnet_mean - swnet_fixed_all_mean))
        bulk_swnet = np.dstack((bulk_swnet, swnet_mean))
        
###############################################################################
# Make DataFrame
###############################################################################

# Make DataFrame
df = pd.DataFrame(list(zip(exp1, exp2, exp3, exp4, exp5)))

df.columns = ['fixed_snow_all', 'observed_albedo', 'fixed_snow_ice', 
              'fixed_ice', 'fixed_snow']

# Divide by area
df = df / np.sum(mask == 1)

# Melt-albedo feedbacks increase SWnet by...
total_albedo = (df['observed_albedo'] - df['fixed_snow_all']) / df['fixed_snow_all']
print(total_albedo.mean())

# Compute radiative forcing in W m-2

# SWnet due to reducing glacier ice albedo 
df['ice_forcing'] = df['observed_albedo'] - df['fixed_ice']

# SWnet due to reducing snow albedo
df['snow_forcing'] = df['observed_albedo'] - df['fixed_snow']

# SWnet due to snowline fluctuations
df['snowline_forcing'] = df['fixed_snow_ice'] - df['fixed_snow_all']

# Save as csv
df.to_csv(path + 'positive-forcing-results.csv')

# Remove first layer
snowline_diff = snowline_diff[:,:,1:]
ice_diff = ice_diff[:,:,1:]
snow_diff = snow_diff[:,:,1:]
bulk_diff = bulk_diff[:,:,1:]
bulk_swnet = bulk_swnet[:,:,1:]

# Save grids as NetCDF
lats, lons = modis['latitude'].values, modis['longitude'].values
    
###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
   
dataset = netCDF4.Dataset(path + 'final-forcing-grids.nc', 'w', format='NETCDF4_CLASSIC')
print('Creating %s' % path + 'final-forcing-grids.nc')
dataset.Title = "Gridded positive radiative forcings"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', snowline_diff.shape[0])
lon_dim = dataset.createDimension('x', snowline_diff.shape[1])
data_dim = dataset.createDimension('z', snowline_diff.shape[2])

    
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
snowline_nc = dataset.createVariable('snowline', np.int8, ('y','x','z'))
snow_nc = dataset.createVariable('snow', np.int8, ('y','x','z'))
ice_nc = dataset.createVariable('ice', np.int8, ('y','x','z'))
bulk_nc = dataset.createVariable('bulk', np.int8, ('y','x','z'))
swnet_nc = dataset.createVariable('swnet', np.int16, ('y','x','z'))

# Write data to layers
Y[:] = lats
X[:] = lons
x[:] = lons[0,:]
y[:] = lats[:,0]
snowline_nc[:] = snowline_diff.astype(np.int8)
snow_nc[:] = snow_diff.astype(np.int8)
ice_nc[:] = ice_diff.astype(np.int8)
bulk_nc[:] = bulk_diff.astype(np.int8)
swnet_nc[:] = bulk_swnet.astype(np.int16)
z[:] = np.arange(1,23)

print('Writing data to %s' % path + 'final-forcing-grids.nc')
    
# Close dataset
dataset.close()




































