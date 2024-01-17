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

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%

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

exp1, exp2, exp3, exp4, exp5 = [], [], [], [], []

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

        # =====================================================================
        # # Fix snow albedo to max physical value
        # # In which case, snow albedo feedback can only be positive
        # fix_snow = np.copy(albedo)
        # fix_snow[classified == 2] = 83
        # 
        # # Compute SWnet
        # swnet_fixed_snow = (1 - (fix_snow / 100)) * merra['swd_allsky']
        # swnet_fixed_snow_mean = np.nanmean(swnet_fixed_snow.values, axis=2)
        #                         
        # # Mask ice sheet
        # swnet_fixed_snow_mean[mask == 0] = np.nan
        # =====================================================================
        
        #######################################################################
        # Fix snow albedo
        #######################################################################
        
        # Fix snow albedo to that observed on June 1
        fix_snow = np.copy(albedo)
        june1 = np.repeat(fix_snow[:,:,0:1], 92, axis=2)
        
        fix_snow[classified == 2] = june1[classified == 2] + j
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
        # Fix both snow and glacier ice albedo
        #######################################################################
        
        # Fix both
        fix_both = np.copy(albedo)
        fix_both[classified == 1] = i
        fix_both[classified == 2] = june1[classified == 2] + j
        
        # Compute SWnet
        swnet_fixed_both = (1 - (fix_both / 100)) * merra['swd_allsky']
        swnet_fixed_both_mean = np.nanmean(swnet_fixed_both.values, axis=2)
                                
        # Mask ice sheet
        swnet_fixed_both_mean[mask == 0] = np.nan
        
        #######################################################################
        # All ice sheet snow
        #######################################################################
        
        # Fix all with snow
        fix_all = np.copy(albedo)
        fix_all[classified > 0] = 83
        
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

#%%
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
print(total_albedo)

#%%

# Save as csv
df.to_csv(path + 'results.csv')

















