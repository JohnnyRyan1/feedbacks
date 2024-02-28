#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute SWnet for negative snowfall feedback experiment

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

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

# Define mask
mask = ismip_1km['GIMP'].values

# Remove exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded_grid_cells.csv')

# Filter
mask[(nan_df['x'].values, nan_df['y'].values)] = False

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%

# Testing
f=0
i=55
lat, lon = 67.0955, -49.951 # KAN-L
#lat, lon = 67.067, -48.836 # KAN-M
#lat, lon = 67.0003,	-47.025 # KAN-U


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

# Classify
classified = np.copy(albedo)
classified[classified <= i] = 1
classified[classified > i] = 2

only_melt = np.copy(albedo)

# First, find the index of the grid point nearest a specific lat/lon.   
abslat = np.abs(modis['latitude'] - lat)
abslon = np.abs(modis['longitude'] - lon)
c = np.maximum(abslon, abslat)

([xloc], [yloc]) = np.where(c == np.min(c))

only_melt = albedo[xloc,yloc,:]

zeros = np.full((only_melt[0].shape), 0)
zeros = zeros[np.newaxis]

only_melt_diff = np.diff(only_melt, axis=0)
only_melt_diff = np.hstack((zeros, only_melt_diff))

only_melt_diff[only_melt_diff < 0] = 0
only_melt_diff[only_melt < 55] = 0
only_melt_cumsum = np.cumsum(only_melt_diff)

melt = np.diff(only_melt, axis=0)
melt = np.hstack((zeros, melt))
melt[melt > 0] = 0

update = only_melt_cumsum + melt

melt_only_new = only_melt - update

# Is it a snow pixel?
if np.mean(classified[xloc,yloc,:]) == 2:
    melt_only_new[melt_only_new <= 55] = 55

if np.mean(classified[xloc, yloc, :]) < 2:
    melt_only_new[melt_only_new <= 30] = 30

# Compute SWnet for observed
swnet = (1 - (only_melt / 100)) * merra['swd_allsky'].values[xloc,yloc,:]
swnet_mean = np.nanmean(swnet, axis=0)

swnet_melt_only = (1 - (melt_only_new / 100)) * merra['swd_allsky'].values[xloc,yloc,:]
swnet_melt_only_mean = np.nanmean(swnet_melt_only, axis=0)


#%%

i = 55

#%%
exp1, exp2 = [], []
diff = np.zeros(mask.shape)

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
            
    # Add min and max snow albedo
    albedo[albedo > 83] = 83
    albedo[albedo <= 30] = 30
    
    # Classify
    classified = np.copy(albedo)
    classified[classified <= i] = 1
    classified[classified > i] = 2
    
    ###########################################################################
    # Perform experiment
    ###########################################################################
    only_melt = np.copy(albedo)
    
    # Add another axis for difference array
    zeros = np.full((albedo[:,:,0].shape), 0)
    zeros = zeros[:,:,np.newaxis]
    
    # Difference along axis
    only_melt_diff = np.diff(only_melt, axis=2)
    only_melt_diff = np.dstack((zeros, only_melt_diff))
    
    # If difference is less than zero (i.e. melt) set to zero
    only_melt_diff[only_melt_diff < 0] = 0
    
    # If albedo is not raised above 55 then set to zero
    only_melt_diff[only_melt < 55] = 0
    only_melt_cumsum = np.cumsum(only_melt_diff, axis=2)
    
    # If difference is more than zero set to zero
    melt = np.diff(only_melt, axis=2)
    melt = np.dstack((zeros, melt))
    melt[melt > 0] = 0
    
    # Add cumulative sum of snowfall to the individual melt events
    update = only_melt_cumsum + melt
    
    # Update array
    only_melt_update = only_melt - update
    
    # Find whether pixel is always snow or not
    snow_classified = (np.nanmean(classified, axis=2) == 2)[:,:,np.newaxis]
    snow = np.repeat(snow_classified, 92, axis=2)
    
    ice_classified = (np.nanmean(classified, axis=2) < 2)[:,:,np.newaxis]
    ice = np.repeat(ice_classified, 92, axis=2)
    
    # Prevent a grid cell in the accumulation from going below 0.56 albedo 
    only_melt_update[(ice == True) & (only_melt_update <= 30)] = 30
    only_melt_update[(snow == True) & (only_melt_update <= 56)] = 56
    
    ###########################################################################
    # Observed albedo
    ###########################################################################

    # Compute SWnet
    swnet = (1 - (albedo / 100)) * merra['swd_allsky']
    swnet_mean = np.nanmean(swnet.values, axis=2)

    # Mask ice sheet
    swnet_mean[mask == 0] = np.nan
    
    ###########################################################################
    # Melt only albedo
    ###########################################################################
    
    # Compute SWnet
    swnet_melt_only = (1 - (only_melt_update / 100)) * merra['swd_allsky']
    swnet_melt_only_mean = np.nanmean(swnet_melt_only.values, axis=2)

    # Mask ice sheet
    swnet_melt_only_mean[mask == 0] = np.nan
    
    # Append to lists
    exp1.append(np.nansum(swnet_mean))
    exp2.append(np.nansum(swnet_melt_only_mean))
    
    # Append to grid
    diff = np.dstack((diff, swnet_melt_only_mean - swnet_mean))

#%%
###############################################################################
# Make DataFrame
###############################################################################

# Make DataFrame
df = pd.DataFrame(list(zip(exp1, exp2)))

df.columns = ['observed_albedo', 'melt_only']

# Divide by area
df = df / np.sum(mask == 1)

# Compute snowfall radiative forcing
df['snowfall_forcing'] = df['melt_only'] - df['observed_albedo']

#%%

# Save as csv
df.to_csv(path + 'negatiave_forcing_results.csv')

#%%

# Save grid as NetCDF














