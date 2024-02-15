#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute empirical relationships between W m-2 and runoff in different elevation bins.

"""

# Import modules
import xarray as xr
import numpy as np
import glob
import matplotlib.pyplot as plt

#%%
# Import data
path = '/Users/jryan4/Dropbox (University of Oregon)/published/snowfall/data/'

# Define MAR files
mar_files = sorted(glob.glob(path + 'mar-v3-12-1/*'))

# Derive MAR surface heights
mar = xr.open_dataset(mar_files[0])
sh = mar['SH'].values

# Define elevations
elevations = np.arange(0, 3400, 200)
#%%

watts_by_elev = []
runoff_by_elev = []
for e in range(len(elevations) - 1):
    elevation_mask = (sh > elevations[e]) & (sh < elevations[e+1])
    
    watts_by_year = []
    runoff_by_year = []
    for m in mar_files:
        mar = xr.open_dataset(m)
        watts = mar['LWD'] + (mar['SWD']*(1-mar['AL2'][:,0,:,:])) + mar['SHF'] + mar['LHF']
        runoff = mar['RU'][:,0,:,:]
                 
        watts_by_year.append(np.nanmean(np.nanmean(watts[152:244], axis=0)[elevation_mask]))
        runoff_by_year.append(np.nanmean(np.nanmean(runoff[152:244], axis=0)[elevation_mask]))
    
    