#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compare MAR vs. downscaled MERRA-2.

"""

# Import modules
import xarray as xr
import glob
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

mar_files = sorted(glob.glob(mar_path + '*.nc'))[2:]
merra = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')


#%%

mar_swd = []
mar_lwd = []

for f in mar_files:
    mar = xr.open_dataset(f)
    mar_lwd_mean = mar['LWD'][151:243, :, :].mean(axis=0).values
    mar_swd_mean = mar['SWD'][151:243, :, :].mean(axis=0).values
    mar_lwd.append(np.nanmean(mar_lwd_mean[mar['MSK'] > 99]))
    mar_swd.append(np.nanmean(mar_swd_mean[mar['MSK'] > 99]))


#%%


merra_lwd = np.mean(merra['lwd_allsky'][:,:,0:20], axis=(0,1))
merra_swd = np.mean(merra['swd_allsky'][:,:,0:20], axis=(0,1))

#%%

plt.scatter(mar_swd, merra_swd)

#%%

plt.scatter(mar_lwd, merra_lwd)


#%%


























