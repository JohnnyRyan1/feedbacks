#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Compute interpolation and gap-filling stats

"""

# Import packages
import xarray as xr
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os

"""
Note that there are 125 values for each grid cell in the interpolated grids but
there are only 92 values for each grid cell in the  gap-filled grids. 

"""

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Define files
files = sorted(glob.glob(path + 'filling-stats/*'))

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define mask
mask = ismip_1km['GIMP'].values

# Import gap-filled results
gap = pd.read_csv(path + 'filling-stats/gap-filled-totals.csv')

#%%

original, cut, interpolated = 0, 0, 0

# Loop over data
for file in files:
    data = pd.read_csv(file)
    original = original + data['original'].sum()
    cut = cut + data['cut'].sum()
    interpolated = interpolated + data['interpolated'].sum()
    
#%%

# Compute total number of daily albedo values for whole ice sheet for 22 years
total_pixels = mask.sum()*125*len(files)

# How many of those grid cells had valid MODIS albedo observations?
print(original/total_pixels)

# How many were removed by the filter?
print(cut/total_pixels)

# How many valid albedo values after interpolation?
print(interpolated/total_pixels)
    
# How many albedo values during summer
summer_pixels = (total_pixels/125) * 92

# Number of albedo values after gap-filling
(summer_pixels - gap['post'].sum())/summer_pixels

#%%
# Total number of exlcuded grid cells
nan_df = pd.read_csv(path + 'excluded_grid_cells.csv')

print(nan_df.shape[0])
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    