#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Interpolate missing data gaps in MODIS albedo data.

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os

#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Import MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo-filled/*.nc'))

# Define destination
dest = path + 'modis-albedo-complete/'

#%%

# Read MODIS data
f = xr.open_dataset(modis_files[0])

# Convert to float
array = f['filled_albedo'].astype(np.float32)

# Replace zeros with NaNs
array = array.where(array != 0)

# Find pixels that contain more than one NaN value
null_pixels = np.sum(np.isnan(array.values), axis=2)
null_pixels = ((ismip_1km['GIMP'] == 1) & (null_pixels > 0))


#%%

# Identify pixel
coords = [np.where(null_pixels)[0][0], np.where(null_pixels)[1][0]]

pixel = f['filled_albedo'][coords[0], coords[1], :]

up = f['filled_albedo'][coords[0]+1, coords[1], :]
down = f['filled_albedo'][coords[0]-1, coords[1], :]
left = f['filled_albedo'][coords[0], coords[1]+1, :]
right = f['filled_albedo'][coords[0], coords[1]-1, :]

up = f['filled_albedo'][coords[0]+1, coords[1]+1, :]
down = f['filled_albedo'][coords[0]-1, coords[1]+1, :]
left = f['filled_albedo'][coords[0]+1, coords[1]-1, :]
right = f['filled_albedo'][coords[0]-1, coords[1]-1, :]











