#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate temperature and cloudiness for ice sheet + regions

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import glob
from scipy import ndimage

#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Covert ISMIP to 50 km grid cells
mask_50km = ndimage.zoom(ismip_1km['GIMP'].values, 0.02)

# Define 3D mask
mask3d = np.repeat(mask_50km[:,:,np.newaxis], 92, axis=2)

# Define regions
regions = xr.open_dataset(path + 'regions-mask.nc')
region_mask = regions['regions'].values
region_mask[region_mask == 0] = np.nan
region_mask = ndimage.zoom(region_mask, 0.02, order=0)

# Define temperature files
t_files = sorted(glob.glob(path + 'merra-t2m-50km/*.nc'))

# Define SWD files
s_files = sorted(glob.glob(path + 'merra-swd-50km/*.nc'))

#%%

t2m_list = []
t2m_region_list = []
swd_allsky_list = []
swd_clrsky_list = []
swd_allsky_region_list = []
swd_clrsky_region_list = []

for j in range(len(t_files)):
    # Import albedo data
    t2m = xr.open_dataset(t_files[j])
    
    # Import albedo data
    swd = xr.open_dataset(s_files[j])
    
    # Apply some preprocessing
    t2m = t2m['t2m'].values.astype(np.float32)
    t2m[mask3d == 0] = np.nan
    t2m[t2m == 0] = np.nan
    
    swd_allsky = swd['swd_allsky'].values.astype(np.float32)
    swd_allsky[mask3d == 0] = np.nan
    swd_allsky[swd_allsky == 0] = np.nan
    
    swd_clrsky = swd['swd_clrsky'].values.astype(np.float32)
    swd_clrsky[mask3d == 0] = np.nan
    swd_clrsky[swd_allsky == 0] = np.nan
    
    # Compute means
    t2m_is = np.nanmean(np.nanmean(t2m, axis=2), axis=(0,1))
    
    t2m_regions = []
    for i in np.arange(1, 9):
        t2m_regions.append(np.nanmean(np.nanmean(t2m, axis=2)[region_mask == i]))
    
    t2m_list.append(t2m_is)
    t2m_region_list.append(t2m_regions)
    
    # Compute means
    swd_allsky_is = np.nanmean(np.nanmean(swd_allsky, axis=2), axis=(0,1))
    swd_clrsky_is = np.nanmean(np.nanmean(swd_clrsky, axis=2), axis=(0,1))
    
    swd_allsky_regions = []
    swd_clrsky_regions = []
    
    for i in np.arange(1, 9):
        swd_allsky_regions.append(np.nanmean(np.nanmean(swd_allsky, axis=2)[region_mask == i]))
        swd_clrsky_regions.append(np.nanmean(np.nanmean(swd_clrsky, axis=2)[region_mask == i]))
    
    swd_allsky_list.append(swd_allsky_is)
    swd_clrsky_list.append(swd_clrsky_is)
    swd_allsky_region_list.append(swd_allsky_regions)
    swd_clrsky_region_list.append(swd_clrsky_regions)
    
# Make DataFrame for ice sheet
df = pd.DataFrame(list(zip(swd_clrsky_list, swd_allsky_list, t2m_list)))
df.columns = ['clrsky', 'allsky', 't2m']
df['cloudiness'] = (df['clrsky'] - df['allsky']) / df['clrsky']
df.to_csv(path + 'ice-sheet-climate.csv')

# Make DataFrame for regions
df_clr = pd.DataFrame(swd_clrsky_region_list)
df_clr.columns = ['1', '2', '3', '4', '5', '6', '7', '8']

df_all = pd.DataFrame(swd_allsky_region_list)
df_all.columns = ['1', '2', '3', '4', '5', '6', '7', '8']

df_cloudiness = (df_clr - df_all) / df_clr

df_t2m = pd.DataFrame(t2m_region_list)
df_t2m.columns = ['1', '2', '3', '4', '5', '6', '7', '8']

# Save as csv
df_cloudiness.to_csv(path + 'regional-cloudiness.csv')
df_t2m.to_csv(path + 'regional-t2m.csv')













