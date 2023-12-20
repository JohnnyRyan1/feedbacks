#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Fill missing data gaps in MODIS albedo data.

"""

# Import modules
import xarray as xr
import numpy as np
import glob

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6.nc')

# Import MODIS files
modis = xr.open_dataset(path + 'modis-albedo/mod10a1-albedo-2002.nc')

#%%

# Find pixel that is in ice mask
modis['albedo'][790,475,:]

# Replace zeros with nan
modis['albedo'] = modis['albedo'].astype(np.float32)
modis['albedo'] = modis['albedo'].where(modis['albedo'] != 0)

# Perform filter (i.e. remove if two standard deviations from mean)
rolling_mean = modis['albedo'][790,475,:].rolling(z=11, min_periods=1, center=True).mean()
rolling_std = modis['albedo'][790,475,:].rolling(z=11, min_periods=1, center=True).std()
rolling_std = rolling_std * 2

# Calculate difference between pixel value and rolling mean
difference = np.abs(modis['albedo'][790,475,:] - rolling_mean)

# Mask values that are more than two standard deviations from the mean
mask = (difference < rolling_std)

# Calculate 11-day rolling median to be used as the timeseries
rolling_median = modis['albedo'][790,475,:].where(mask == True).rolling(z=11, min_periods=3, center=True).median()

# Linearly interpolate between values
linear_interp = rolling_median.interpolate_na(dim="z", method="linear")

#%%
# Add to original Dataset
modis['albedo_gf'] = linear_interp

# Classify
modis['ice'] = ((modis['albedo_gf'] <= 55) & (modis['albedo_gf'] != 0)).astype(np.int8)


# Perform filter (i.e. remove if two standard deviations from mean)
rolling_mean = albedo.rolling(z=11, min_periods=1, center=True).mean()
rolling_std = albedo.rolling(z=11, min_periods=1, center=True).std()
rolling_std = rolling_std * 2

# Calculate difference between pixel value and rolling mean
difference = np.abs(albedo - rolling_mean)

# Mask values that are more than two standard deviations from the mean
mask = (difference < rolling_std)

# Calculate 11-day rolling median to be used as the timeseries
rolling_median = albedo.where(mask == True).rolling(z=11, min_periods=3, center=True).median()








