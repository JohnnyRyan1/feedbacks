#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute SWnet for clouds experiments

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import netCDF4
import itertools

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define MERRA
merra = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')

# Define MERRA climatology
modis_climatology = xr.open_dataset(path + 'modis-climatology.nc')
albedo = modis_climatology['albedo_climatology'].values.astype('float')
swd_lwd = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-melt-coeffs.csv')

#%%

# Define mask
mask = ismip_1km['GIMP'].values

# Define ice threshold
ice_threshold = [53, 55, 57]

# Define albedo uncertainty
snow_threshold = [84]

#%%

# Remove exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded-grid-cells.csv')

# Filter
mask[(nan_df['x'].values, nan_df['y'].values)] = False

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%
# Some preprocessing
albedo[mask3d == 0] = np.nan
albedo[albedo == 0] = np.nan

exp1, exp2 = [], []
sw_bulk_cloud = np.zeros(mask.shape)
sw_bulk_no_cloud = np.zeros(mask.shape)

for i, j in itertools.product(ice_threshold, snow_threshold):

    # Add max snow albedo
    albedo[albedo > j] = j
    albedo[albedo <= 30] = 30
    
    # Classify
    classified = np.copy(albedo)
    classified[classified <= i] = 1
    classified[classified > i] = 2
    
    observed = np.copy(albedo)
    observed = observed
    
    # Compute mean climatology
    observed_mean = np.nanmean(observed, axis=2)
    
    for f in range(merra['z'].shape[0]):
    
        print('Processing... %s' % f)
      
        # Compute SWnet with observed clouds
        swnet_mean_cloud = (1 - (observed_mean / 100)) * merra['swd_allsky'][:,:,f].values
    
        # Mask ice sheet
        swnet_mean_cloud[mask == 0] = np.nan
        
        # Compute SWnet with no clouds
        swnet_mean_nocloud = (1 - (observed_mean / 100)) * merra['swd_clrsky'][:,:,f].values
    
        # Mask ice sheet
        swnet_mean_nocloud[mask == 0] = np.nan
        
        # Append to lists
        exp1.append(np.nansum(swnet_mean_cloud))
        exp2.append(np.nansum(swnet_mean_nocloud))
       
        # Append to grids
        sw_bulk_cloud = np.dstack((sw_bulk_cloud, swnet_mean_cloud))
        sw_bulk_no_cloud = np.dstack((sw_bulk_no_cloud, swnet_mean_nocloud))
        
# Remove first layer
sw_bulk_cloud = sw_bulk_cloud[:,:,1:]
sw_bulk_no_cloud = sw_bulk_no_cloud[:,:,1:]

#%%


# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

clouds_sw_low, clouds_sw_med, clouds_sw_high, clouds_lw = [], [], [], []
area = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    
    cre_sw_low = np.nanmean((sw_bulk_no_cloud[:,:,0:22] - sw_bulk_cloud[:,:,0:22]), axis=2)
    cre_sw_med = np.nanmean((sw_bulk_no_cloud[:,:,22:44] - sw_bulk_cloud[:,:,22:44]), axis=2)
    cre_sw_high = np.nanmean((sw_bulk_no_cloud[:,:,44:] - sw_bulk_cloud[:,:,44:]), axis=2)
    cre_lw_mean = np.nanmean((swd_lwd['lwd_allsky'].values - swd_lwd['lwd_clrsky'].values), axis=2)
    
    clouds_sw_low.append(np.nanmean(cre_sw_low[elevation_mask]))
    clouds_sw_med.append(np.nanmean(cre_sw_med[elevation_mask]))
    clouds_sw_high.append(np.nanmean(cre_sw_high[elevation_mask]))

    clouds_lw.append(np.nanmean(cre_lw_mean[elevation_mask]))
                         
area = np.array(area)

#%%

# Compute melt factor
coeffs['factor'] = coeffs['melt'] / coeffs['watts']

coeffs['clouds_sw_low'] = coeffs['factor'] * clouds_sw_low
coeffs['clouds_sw_med'] = coeffs['factor'] * clouds_sw_med
coeffs['clouds_sw_high'] = coeffs['factor'] * clouds_sw_high

coeffs['clouds_lw'] = coeffs['factor'] * clouds_lw


coeffs['cloud_sw_gt_low'] = coeffs['clouds_sw_low'] / 1e+06 * area * 92
coeffs['cloud_sw_gt_med'] = coeffs['clouds_sw_med'] / 1e+06 * area * 92
coeffs['cloud_sw_gt_high'] = coeffs['clouds_sw_high'] / 1e+06 * area * 92

coeffs['cloud_lw_gt'] = coeffs['clouds_lw'] / 1e+06 * area * 92



#%%

# Reduction of meltwater melt due to clouds
print(coeffs['cloud_lw_gt'].sum() - coeffs['cloud_sw_gt_low'].sum())
print(coeffs['cloud_lw_gt'].sum() - coeffs['cloud_sw_gt_high'].sum())

print(np.sqrt(coeffs['cloud_lw_gt_std'].sum()**2 + coeffs['cloud_sw_gt_std'].sum()**2))




























