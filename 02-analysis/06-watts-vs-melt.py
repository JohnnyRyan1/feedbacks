#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute empirical relationships between W m-2 and melt in different elevation bins.

"""

# Import modules
import xarray as xr
import numpy as np
import glob
from scipy import stats
import pandas as pd

#%%

# Define user
user = 'jryan4'

# Import data
alt_path = '/Users/jryan4/Dropbox (University of Oregon)/published/snowfall/data/'
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Define MAR files
mar_files = sorted(glob.glob(alt_path + 'mar-v3-12-1/*'))

# Derive MAR surface heights
mar = xr.open_dataset(mar_files[0])
sh = mar['SH'].values
msk = mar['MSK'].values

# Define elevations
elevations = np.arange(0, 3600, 200)
#%%

watts_by_elev = []
melt_by_elev = []
runoff_by_elev = []
area_by_elev = []
for e in range(len(elevations) - 1):
    print(e)
    elevation_mask = (sh > elevations[e]) & (sh < elevations[e+1]) & (msk > 90)
    
    watts_by_year = []
    melt_by_year = []
    runoff_by_year = []
    for m in mar_files:
        mar = xr.open_dataset(m)
        watts = mar['LWD'] + (mar['SWD']*(1-mar['AL2'][:,0,:,:])) + mar['SHF'] + mar['LHF']
        melt = mar['ME'][:,0,:,:]     
        runoff = mar['RU'][:,0,:,:]      
                 
        watts_by_year.append(np.nanmean(np.nanmean(watts[152:244, :, :], axis=0)[elevation_mask]))
        melt_by_year.append(np.nanmean(np.nanmean(melt[152:244, :, :], axis=0)[elevation_mask]))
        runoff_by_year.append(np.nanmean(np.nanmean(runoff[152:244, :, :], axis=0)[elevation_mask]))

    watts_by_elev.append(watts_by_year)
    runoff_by_elev.append(runoff_by_year)
    melt_by_elev.append(melt_by_year)
    area_by_elev.append(np.sum(elevation_mask)*10)
        
#%%

# Save coefficients as csv
coeffs = []
for i in range(len(elevations) - 1):
    coeffs.append(stats.linregress(watts_by_elev[i], melt_by_elev[i]))

watts_df = pd.DataFrame(watts_by_elev)
melt_df = pd.DataFrame(melt_by_elev)

watts_by_elev = np.array(watts_by_elev)
melt_by_elev = np.array(melt_by_elev)
df1 = pd.DataFrame(melt_by_elev/watts_by_elev)
df2 = pd.DataFrame(melt_by_elev)
df3 = pd.DataFrame(watts_by_elev)
df4 = pd.DataFrame(runoff_by_elev)

df = pd.DataFrame(coeffs)
df['elev'] = elevations[:-1]
df['watts'] = watts_df.mean(axis=1)
df['melt'] = melt_df.mean(axis=1)

# Save as csv
df.to_csv(path + 'watts-to-melt-coeffs.csv')
df1.to_csv(path + 'melt-factors.csv', index=False)
df2.to_csv(path + 'melt.csv', index=False)
df3.to_csv(path + 'watts.csv', index=False)
df4.to_csv(path + 'runoff.csv', index=False)


#%%


# Temperature vs. ice sheet melt

temp_by_year = []
melt_by_year = []
for m in mar_files:
    mar = xr.open_dataset(m)
    melt = mar['ME'][:,0,:,:]
    temp = mar['TT'][:,0,:,:]
    
    melt_by_year.append(np.nansum(np.nansum(melt[152:244, :, :], axis=0)[(msk > 90)]))
    temp_by_year.append(np.nanmean(np.nanmean(temp[152:244, :, :], axis=0)[(msk > 90)]))


melt_by_year = np.array(melt_by_year)
temp_by_year = np.array(temp_by_year)

# Convert from mm to Gt
melt_gt = melt_by_year / 1e+6 * 10 * 10
temp_anom = temp_by_year - np.mean(temp_by_year)

df = pd.DataFrame(list(zip(temp_anom, melt_gt)))
df.columns = ['temp', 'melt']
df.to_csv(path + 'mar-vs-temp.csv')

#%%

    
    
    
    
    
    
    
    
    
    