#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Raise in albedo due to cloud correction.

"""

# Import modules
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values


# Import data
post = xr.open_dataset(path + 'modis-albedo-final/mod10a1-albedo-2020.nc')
pre = xr.open_dataset(path + 'modis-albedo-complete/mod10a1-albedo-2020.nc')

#%%

mean_pre = np.nanmean(np.nanmean(pre['albedo'], axis=2)[mask])

mean_post = np.nanmean(np.nanmean(post['albedo'], axis=2)[mask])


diff = np.nanmean(pre['albedo'], axis=2)[mask] - np.nanmean(post['albedo'], axis=2)[mask]





