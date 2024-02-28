#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute melt-albedo feedback at every grid cell

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
from scipy import stats
import pandas as pd

"""
NOTE: Can only compute feedback strength for snow over pixels that were always snow during the 
study period.
"""
#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-runoff-coeffs.csv')

#%%

# Import data
data = xr.open_dataset(path + 'final-forcing-grids.nc')
t2m = xr.open_dataset(path + 'final-temp-grids.nc')
swd = xr.open_dataset(path + 'final-swd-grids.nc')
snow_mask = xr.open_dataset(path + 'snow-mask.nc')
cloud = xr.open_dataset(path +'final-cloud-forcing-grids.nc')

# Compute anomalies
data['swnet_anom'] = data['swnet'] - data['swnet'].mean(dim='z')
t2m['t2m_anom'] = t2m['t2m'] - t2m['t2m'].mean(dim='z')
swd['cre'] = cloud['swnet_no_cloud'] - cloud['swnet_cloud']
swd['cre_anom'] = swd['cre'] - swd['cre'].mean(dim='z')

#%%
# Define some functions
def new_linregress(x, y):
    # Wrapper around scipy linregress to use in apply_ufunc
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return np.array([slope, intercept, r_value, p_value, std_err])


#%%
###############################################################################
# Compute linear relationship between radiative forcing and air temperature for every pixel
###############################################################################
xr_stats_temp = xr.apply_ufunc(new_linregress, t2m['t2m_anom'], 
                          data['swnet_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_cloud = xr.apply_ufunc(new_linregress, t2m['t2m_anom'], 
                          swd['cre_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})


xr_stats_snowline = xr.apply_ufunc(new_linregress, t2m['t2m_anom'], 
                          data['snowline'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_snow = xr.apply_ufunc(new_linregress, t2m['t2m_anom'], 
                          data['snow'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_ice = xr.apply_ufunc(new_linregress, t2m['t2m_anom'], 
                          data['ice'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

#%%

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'feedback-stats.nc', 
                          'w', format='NETCDF4_CLASSIC')
print('Creating... %s' % path + 'feedback-stats.nc')
dataset.Title = "Slopes and intercepts for temperature vs. albedo relationship"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C., (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', xr_stats_temp.shape[0])
lon_dim = dataset.createDimension('x', xr_stats_temp.shape[1])
data_dim =  dataset.createDimension('z', xr_stats_temp.shape[2])

# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
xr_stats_temp_nc = dataset.createVariable('bulk_vs_temp', np.float32, ('y','x','z'))
xr_stats_cloud_nc = dataset.createVariable('bulk_vs_cloud', np.float32, ('y','x','z'))

xr_stats_snowline_nc = dataset.createVariable('snowline', np.float32, ('y','x','z'))
xr_stats_snow_nc = dataset.createVariable('snow', np.float32, ('y','x','z'))
xr_stats_ice_nc = dataset.createVariable('ice', np.float32, ('y','x','z'))

# Write data to layers
Y[:] = ismip_1km['lat'].values
X[:] = ismip_1km['lon'].values
xr_stats_temp_nc[:] = xr_stats_temp.values
xr_stats_cloud_nc[:] = xr_stats_cloud.values
xr_stats_snowline_nc[:] = xr_stats_snowline.values
xr_stats_snow_nc[:] = xr_stats_snow.values
xr_stats_ice_nc[:] = xr_stats_ice.values

print('Writing data to %s' % path + 'feedback-stats.nc')
    
# Close dataset
dataset.close()




















