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
modis_files = sorted(glob.glob(path + 'modis-albedo-final/*.nc'))

# Define MERRA files
merra = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')

# Compute SW climatology
merra_climatology = merra['swd_allsky'].mean(axis=2).values
merra_climatology[merra_climatology == 0] = np.nan

#%%

# Define ice threshold
ice_threshold = [53, 55, 57]

# Define albedo uncertainty
snow_threshold = [0]

# Define mask
mask = ismip_1km['GIMP'].values

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)


#%%

# Remove exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded-grid-cells.csv')

# Filter
mask[(nan_df['x'].values, nan_df['y'].values)] = False

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 92, axis=2)

#%%

exp1, exp2, exp3, exp4, exp5 = [], [], [], [], []
ice_diff, snow_diff, snowline_diff = np.zeros(mask.shape), np.zeros(mask.shape), np.zeros(mask.shape)
bulk_diff, bulk_swnet = np.zeros(mask.shape), np.zeros(mask.shape)

for i, j in itertools.product(ice_threshold, snow_threshold):

    for f in range(len(modis_files)):

        print('Processing... %s' % f)

        # Import albedo data
        modis = xr.open_dataset(modis_files[f])

        # Some preprocessing
        albedo = modis['albedo'].values.astype(np.float32)
        albedo[mask3d == 0] = np.nan
        albedo[albedo == 0] = np.nan
                
        # Add max snow albedo
        albedo[albedo > 84] = 84
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
        swnet_mean = (1 - (observed.mean(axis=2) / 100)) * merra_climatology

        # Mask ice sheet
        swnet_mean[mask == 0] = np.nan
               
        #######################################################################
        # Fix snow albedo
        #######################################################################
        
        # Fix snow albedo to fresh snow
        fix_snow = np.copy(albedo)       
        
        fix_snow[classified == 2] = 84 + j
        fix_snow[classified == 1] = fix_snow[classified == 1] + j
        
        # Compute SWnet
        swnet_fixed_snow_mean = (1 - (fix_snow.mean(axis=2) / 100)) * merra_climatology
                               
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
        swnet_fixed_ice_mean = (1 - (fix_ice.mean(axis=2) / 100)) * merra_climatology
                                
        # Mask ice sheet
        swnet_fixed_ice_mean[mask == 0] = np.nan
        
        #######################################################################
        # Fix snow albedo to that 84 and glacier ice albedo to 55
        #######################################################################
        
        # Fix both
        fix_both = np.copy(albedo)
        fix_both[classified == 1] = i
        fix_both[classified == 2] = 84 + j
        
        # Compute SWnet
        swnet_fixed_both_mean = (1 - (fix_both.mean(axis=2) / 100)) * merra_climatology
                                
        # Mask ice sheet
        swnet_fixed_both_mean[mask == 0] = np.nan
        
        #######################################################################
        # Fix all the ice sheet to 84
        #######################################################################
        
        # Fix all with snow and glacier ice albedo observed on June 1
        fix_all = np.copy(albedo)
        fix_all[classified > 0] = 84 + j
        
        # Compute SWnet
        swnet_fixed_all_mean = (1 - (fix_all.mean(axis=2) / 100)) * merra_climatology
        
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
        
# Remove first layer
snowline_diff = snowline_diff[:,:,1:]
ice_diff = ice_diff[:,:,1:]
snow_diff = snow_diff[:,:,1:]
bulk_diff = bulk_diff[:,:,1:]
bulk_swnet = bulk_swnet[:,:,1:]

   
#%%

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

ice_low, ice_med, ice_high = [] , [], []
snowline_low, snowline_med, snowline_high = [], [], []
snow_low, snow_med, snow_high = [], [], []
swnet_low, swnet_med, swnet_high = [], [], []
area = []

for e in range(len(elevations) - 1):
    print(e)
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    
    ice_low.append(np.nanmean(np.nanmean(np.nanmean(ice_diff[:,:,0:22], axis=2)[elevation_mask])))
    ice_med.append(np.nanmean(np.nanmean(np.nanmean(ice_diff[:,:,22:44], axis=2)[elevation_mask])))
    ice_high.append(np.nanmean(np.nanmean(np.nanmean(ice_diff[:,:,44:], axis=2)[elevation_mask])))

    snowline_low.append(np.nanmean(np.nanmean(np.nanmean(snowline_diff[:,:,0:22], axis=2)[elevation_mask])))
    snowline_med.append(np.nanmean(np.nanmean(np.nanmean(snowline_diff[:,:,22:44], axis=2)[elevation_mask])))
    snowline_high.append(np.nanmean(np.nanmean(np.nanmean(snowline_diff[:,:,44:], axis=2)[elevation_mask])))
    
    snow_low.append(np.nanmean(np.nanmean(np.nanmean(snow_diff[:,:,0:22], axis=2)[elevation_mask])))
    snow_med.append(np.nanmean(np.nanmean(np.nanmean(snow_diff[:,:,22:44], axis=2)[elevation_mask])))
    snow_high.append(np.nanmean(np.nanmean(np.nanmean(snow_diff[:,:,44:], axis=2)[elevation_mask])))
    
    swnet_low.append(np.nanmean(np.nanmean(np.nanmean(bulk_swnet[:,:,0:22], axis=2)[elevation_mask])))
    swnet_med.append(np.nanmean(np.nanmean(np.nanmean(bulk_swnet[:,:,22:44], axis=2)[elevation_mask])))
    swnet_high.append(np.nanmean(np.nanmean(np.nanmean(bulk_swnet[:,:,44:], axis=2)[elevation_mask])))


area = np.array(area)

#%%

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-melt-coeffs.csv')
factor = pd.read_csv(path + 'melt-factors.csv')

# Compute melt factor
coeffs['factor'] = coeffs['melt'] / coeffs['watts']

coeffs['ice_low'] = coeffs['factor'] * ice_low
coeffs['ice_med'] = coeffs['factor'] * ice_med
coeffs['ice_high'] = coeffs['factor'] * ice_high

coeffs['snowline_low'] = coeffs['factor'] * snowline_low
coeffs['snowline_med'] = coeffs['factor'] * snowline_med
coeffs['snowline_high'] = coeffs['factor'] * snowline_high

coeffs['snow_low'] = coeffs['factor'] * snow_low
coeffs['snow_med'] = coeffs['factor'] * snow_med
coeffs['snow_high'] = coeffs['factor'] * snow_high

# Convert to Gt (mm to km, multiply by area, multiply by number of days)
coeffs['melt_gt'] = coeffs['melt']/1e+06*area*92

coeffs['ice_low_gt'] = coeffs['ice_low']/1e+06*area*92
coeffs['snowline_low_gt'] = coeffs['snowline_low']/1e+06*area*92
coeffs['snow_low_gt'] = coeffs['snow_low']/1e+06*area*92

coeffs['ice_med_gt'] = coeffs['ice_med']/1e+06*area*92
coeffs['snowline_med_gt'] = coeffs['snowline_med']/1e+06*area*92
coeffs['snow_med_gt'] = coeffs['snow_med']/1e+06*area*92

coeffs['ice_high_gt'] = coeffs['ice_high']/1e+06*area*92
coeffs['snowline_high_gt'] = coeffs['snowline_high']/1e+06*area*92
coeffs['snow_high_gt'] = coeffs['snow_high']/1e+06*area*92


#%%
"""
P1
"""

# Amount of meltwater melt due to 0.02 change in glacier ice albedo threshold
print(np.sum(coeffs['snowline_med_gt']), '+/-', np.sum(coeffs['snowline_low_gt']),
                                                np.sum(coeffs['snowline_high_gt']))

print(np.sum(coeffs['snow_med_gt']), '+/-', np.sum(coeffs['snow_low_gt']),
                                                np.sum(coeffs['snow_high_gt']))

print(np.sum(coeffs['ice_med_gt']), '+/-', np.sum(coeffs['ice_low_gt']),
                                                np.sum(coeffs['ice_high_gt']))

# Change per 0.01 in albedo threshold as a percentage
print((np.sum(coeffs['ice_med_gt']) - np.sum(coeffs['ice_low_gt'])) / np.sum(coeffs['ice_low_gt'])/2)
print((np.sum(coeffs['snowline_med_gt']) - np.sum(coeffs['snowline_low_gt'])) / np.sum(coeffs['snowline_low_gt'])/2)
print((np.sum(coeffs['snow_med_gt']) - np.sum(coeffs['snow_low_gt'])) / np.sum(coeffs['snow_low_gt'])/2)

# Bulk forcing is the same
all_low_gt = np.sum(coeffs['ice_low_gt']) + np.sum(coeffs['snowline_low_gt']) + np.sum(coeffs['snow_low_gt'])
all_med_gt = np.sum(coeffs['ice_med_gt']) + np.sum(coeffs['snowline_med_gt']) + np.sum(coeffs['snow_med_gt'])
all_high_gt = np.sum(coeffs['ice_high_gt']) + np.sum(coeffs['snowline_high_gt']) + np.sum(coeffs['snow_high_gt'])
               
# Contribution of glacier ice albedo to bulk radiative forcing
print(np.sum(coeffs['ice_high_gt']) / 
      (np.sum(coeffs['ice_high_gt']) + np.sum(coeffs['snowline_high_gt']) + np.sum(coeffs['snow_high_gt'])))

print(np.sum(coeffs['ice_med_gt']) / 
      (np.sum(coeffs['ice_med_gt']) + np.sum(coeffs['snowline_med_gt']) + np.sum(coeffs['snow_med_gt'])))

print(np.sum(coeffs['ice_low_gt']) / 
      (np.sum(coeffs['ice_low_gt']) + np.sum(coeffs['snowline_low_gt']) + np.sum(coeffs['snow_low_gt'])))


# Contribution of snowline albedo to bulk radiative forcing
print(np.sum(coeffs['snowline_high_gt']) / 
      (np.sum(coeffs['ice_high_gt']) + np.sum(coeffs['snowline_high_gt']) + np.sum(coeffs['snow_high_gt'])))

print(np.sum(coeffs['snowline_med_gt']) / 
      (np.sum(coeffs['ice_med_gt']) + np.sum(coeffs['snowline_med_gt']) + np.sum(coeffs['snow_med_gt'])))

print(np.sum(coeffs['snowline_low_gt']) / 
      (np.sum(coeffs['ice_low_gt']) + np.sum(coeffs['snowline_low_gt']) + np.sum(coeffs['snow_low_gt'])))


# Contribution of snow albedo to bulk radiative forcing
print(np.sum(coeffs['snow_high_gt']) / 
      (np.sum(coeffs['ice_high_gt']) + np.sum(coeffs['snowline_high_gt']) + np.sum(coeffs['snow_high_gt'])))

print(np.sum(coeffs['snow_med_gt']) / 
      (np.sum(coeffs['ice_med_gt']) + np.sum(coeffs['snowline_med_gt']) + np.sum(coeffs['snow_med_gt'])))

print(np.sum(coeffs['snow_low_gt']) / 
      (np.sum(coeffs['ice_low_gt']) + np.sum(coeffs['snowline_low_gt']) + np.sum(coeffs['snow_low_gt'])))



















