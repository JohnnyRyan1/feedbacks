#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare MAR vs. downscaled MERRA-2.

"""

# Import modules
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import rioxarray

# Define user
user = 'johnnyryan'

# Define path
mar_path = '/Users/' + user + '/Dropbox (University of Oregon)/published/snowfall/data/mar-v3-12-1/'
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

#%%

mar = xr.open_dataset(mar_path + 'MARv3.12.1-10km-daily-ERA5-2019.nc')
merra = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')


#%%

mar_mean = mar['LWD'][151:243, :, :].mean(axis=0)

merra_mean_lwd = merra['lwd_allsky'][:,:,17]
merra_mean_swd = merra['swd_allsky'][:,:,17]

mask = mar['MSK'] > 99

#%%


























