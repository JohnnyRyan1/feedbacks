#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Calculate temperature and cloudiness for ice sheet + regions

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import glob
from scipy import ndimage

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Get ISMIP6 lat lons
lon_1km = ismip_1km.variables['lon'][:]
lat_1km = ismip_1km.variables['lat'][:]

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

#%%

"""
Convert daily air temperatures to summer air temperatures

"""


# Define temperature files
files = sorted(glob.glob(path + 'merra-t2m-resample/*'))

t2m = np.zeros((2881,1681))
for f in files:
    t2m_data = xr.open_dataset(f)
    t2m_is = np.nanmean(t2m_data['t2m'], axis=2)
    t2m = np.dstack((t2m, t2m_is))

t2m = t2m[:,:,1:]

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'final-temp-grids.nc', 
                          'w', format='NETCDF4_CLASSIC')
print('Creating %s' %path + 'final-temp-grids.nc')
dataset.Title = "Summer air temperature from MERRA-2"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', t2m.shape[0])
lon_dim = dataset.createDimension('x', t2m.shape[1])
data_dim = dataset.createDimension('z', t2m.shape[2])
    
# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))

y = dataset.createVariable('y', np.float32, ('y'))
x = dataset.createVariable('x', np.float32, ('x'))
z = dataset.createVariable('z', np.float32, ('z'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
t2m_nc = dataset.createVariable('t2m', np.float32, ('y','x','z'))

# Write data to layers
Y[:] = lat_1km
X[:] = lon_1km
x[:] = lon_1km[0,:]
y[:] = lat_1km[:,0]
t2m_nc[:] = t2m.astype(np.float32)
z[:] = np.arange(1,23)

print('Writing data to %s' %path + 'final-temp-grids.nc')
    
# Close dataset
dataset.close()

#%%

"""
Convert SWD to mean summer values

"""

# Define temperature files
files = sorted(glob.glob(path + 'merra-swd-resample/*'))

swd_all = np.zeros((2881,1681))
swd_clr = np.zeros((2881,1681))

for f in files:
    swd_data = xr.open_dataset(f)
    swd__all_is = np.nanmean(swd_data['swd_allsky'], axis=2)
    swd_all = np.dstack((swd_all, swd__all_is))
    swd__clr_is = np.nanmean(swd_data['swd_clrsky'], axis=2)
    swd_clr = np.dstack((swd_clr, swd__clr_is))

swd_all = swd_all[:,:,1:]
swd_clr = swd_clr[:,:,1:]

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'final-swd-grids.nc', 
                          'w', format='NETCDF4_CLASSIC')
print('Creating %s' %path + 'final-swd-grids.nc')
dataset.Title = "Summer SWD from MERRA-2"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', swd_all.shape[0])
lon_dim = dataset.createDimension('x', swd_all.shape[1])
data_dim = dataset.createDimension('z', swd_all.shape[2])
    
# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))

y = dataset.createVariable('y', np.float32, ('y'))
x = dataset.createVariable('x', np.float32, ('x'))
z = dataset.createVariable('z', np.float32, ('z'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
swd_all_nc = dataset.createVariable('swd_all', np.float32, ('y','x','z'))
swd_clr_nc = dataset.createVariable('swd_clr', np.float32, ('y','x','z'))

# Write data to layers
Y[:] = lat_1km
X[:] = lon_1km
x[:] = lon_1km[0,:]
y[:] = lat_1km[:,0]
swd_all_nc[:] = swd_all.astype(np.float32)
swd_clr_nc[:] = swd_clr.astype(np.float32)
z[:] = np.arange(1,23)

print('Writing data to %s' %path + 'final-swd-grids.nc')
    
# Close dataset
dataset.close()




