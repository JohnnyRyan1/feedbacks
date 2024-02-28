#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute SWnet for clouds experiments

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

# Define MERRA files
merra_files = sorted(glob.glob(path + 'merra-swd-resample/*.nc'))

# Define MERRA climatology
m = xr.open_dataset(path + 'modis-climatology.nc')

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
i = 55

# Define albedo uncertainty
j = 0

#%%

# Remove exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded_grid_cells.csv')

# Filter
mask[(nan_df['x'].values, nan_df['y'].values)] = False

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)


#%%

exp1, exp2 = [], []
bulk_cloud = np.zeros(mask.shape)
bulk_no_cloud = np.zeros(mask.shape)

for f in range(len(merra_files)):

    print('Processing... %s' % merra_files[f])

    # Import SW data
    merra = xr.open_dataset(merra_files[f])

    # Some preprocessing
    albedo = m['albedo_climatology'].values.astype(np.float32)
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
   
    # Compute SWnet with observed clouds
    swnet_cloud = (1 - (observed / 100)) * merra['swd_allsky']
    swnet_mean_cloud = np.nanmean(swnet_cloud.values, axis=2)

    # Mask ice sheet
    swnet_mean_cloud[mask == 0] = np.nan
    
    # Compute SWnet with no clouds
    swnet_nocloud = (1 - (observed / 100)) * merra['swd_clrsky']
    swnet_mean_nocloud = np.nanmean(swnet_nocloud.values, axis=2)

    # Mask ice sheet
    swnet_mean_nocloud[mask == 0] = np.nan
    
    # Append to lists
    exp1.append(np.nansum(swnet_mean_cloud))
    exp2.append(np.nansum(swnet_mean_nocloud))
   
    # Append to grids
    bulk_cloud = np.dstack((bulk_cloud, swnet_mean_cloud))
    bulk_no_cloud = np.dstack((bulk_no_cloud, swnet_mean_nocloud))
        
###############################################################################
# Make DataFrame
###############################################################################

# Make DataFrame
df = pd.DataFrame(list(zip(exp1, exp2)))

df.columns = ['cloud', 'no_cloud']

# Divide by area
df = df / np.sum(mask == 1)

# Save as csv
df.to_csv(path + 'cloud-forcing-results.csv')

# Remove first layer
bulk_cloud = bulk_cloud[:,:,1:]
bulk_no_cloud = bulk_no_cloud[:,:,1:]

# Save grids as NetCDF
lats, lons = m['latitude'].values, m['longitude'].values
    
###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
   
dataset = netCDF4.Dataset(path + 'final-cloud-forcing-grids.nc', 'w', format='NETCDF4_CLASSIC')
print('Creating %s' % path + 'final-cloud-forcing-grids.nc')
dataset.Title = "Gridded radiative forcings for clouds"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', bulk_cloud.shape[0])
lon_dim = dataset.createDimension('x', bulk_cloud.shape[1])
data_dim = dataset.createDimension('z', bulk_cloud.shape[2])

    
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
cloud_nc = dataset.createVariable('swnet_cloud', np.int16, ('y','x','z'))
no_cloud_nc = dataset.createVariable('swnet_no_cloud', np.int16, ('y','x','z'))

# Write data to layers
Y[:] = lats
X[:] = lons
x[:] = lons[0,:]
y[:] = lats[:,0]
cloud_nc[:] = bulk_cloud.astype(np.int16)
no_cloud_nc[:] = bulk_no_cloud.astype(np.int16)
z[:] = np.arange(1,23)

print('Writing data to %s' % path + 'final-cloud-forcing-grids.nc')
    
# Close dataset
dataset.close()




































