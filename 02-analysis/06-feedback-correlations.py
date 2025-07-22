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

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
alt_path = '/Users/' + user + '/Dropbox (University of Oregon)/published/clouds/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-runoff-coeffs.csv')

#%%

# Import data
surface = xr.open_dataset(path + 'final-surface-forcing-grids.nc')
cloud = xr.open_dataset(path + 'final-cloud-forcing-grids.nc')
t2m_lwd_swd = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')

# Compute anomalies
surface['swnet_anom'] = surface['swnet'] - surface['swnet'].mean(dim='z')

cloud['swnet_diff'] = cloud['swnet_cloud'] - cloud['swnet_no_cloud']
cloud['swnet_anom'] = cloud['swnet_diff'] - cloud['swnet_diff'].mean(dim='z')

t2m_lwd_swd['t2m_anom'] = t2m_lwd_swd['t2m'] - t2m_lwd_swd['t2m'].mean(dim='z')

t2m_lwd_swd['lwd'] = t2m_lwd_swd['lwd_allsky'] - t2m_lwd_swd['lwd_clrsky']
t2m_lwd_swd['lwd_anom'] = t2m_lwd_swd['lwd'] - t2m_lwd_swd['lwd'].mean(dim='z')

t2m_lwd_swd['net'] = t2m_lwd_swd['lwd'] + cloud['swnet_diff']
t2m_lwd_swd['net_anom'] = t2m_lwd_swd['net'] - t2m_lwd_swd['net'].mean(dim='z')

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
xr_stats_surf = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          surface['swnet_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_sw = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          cloud['swnet_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})


xr_stats_lw = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          t2m_lwd_swd['lwd_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})


xr_stats_net = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          t2m_lwd_swd['net_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

#%%
xr_stats_snowline = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          surface['snowline'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_snow = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          surface['snow'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_ice = xr.apply_ufunc(new_linregress, t2m_lwd_swd['t2m_anom'], 
                          surface['ice'],
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
lat_dim = dataset.createDimension('y', xr_stats_surf.shape[0])
lon_dim = dataset.createDimension('x', xr_stats_surf.shape[1])
data_dim =  dataset.createDimension('z', xr_stats_surf.shape[2])

# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
xr_stats_surf_nc = dataset.createVariable('surface_feedback', np.float32, ('y','x','z'))
xr_stats_lw_nc = dataset.createVariable('cloud_lw_feedback', np.float32, ('y','x','z'))
xr_stats_sw_nc = dataset.createVariable('cloud_sw_feedback', np.float32, ('y','x','z'))
xr_stats_net_nc = dataset.createVariable('cloud_net_feedback', np.float32, ('y','x','z'))

xr_stats_snowline_nc = dataset.createVariable('snowline_feedback', np.float32, ('y','x','z'))
xr_stats_snow_nc = dataset.createVariable('snow_feedback', np.float32, ('y','x','z'))
xr_stats_ice_nc = dataset.createVariable('ice_feedback', np.float32, ('y','x','z'))

# Write data to layers
Y[:] = ismip_1km['lat'].values
X[:] = ismip_1km['lon'].values
xr_stats_surf_nc[:] = xr_stats_surf.values
xr_stats_lw_nc[:] = xr_stats_lw.values
xr_stats_sw_nc[:] = xr_stats_sw.values
xr_stats_net_nc[:] = xr_stats_net.values
xr_stats_snowline_nc[:] = xr_stats_snowline.values
xr_stats_snow_nc[:] = xr_stats_snow.values
xr_stats_ice_nc[:] = xr_stats_ice.values

print('Writing data to %s' % path + 'feedback-stats.nc')
    
# Close dataset
dataset.close()

#%%

# Checks
xr_stats_sw_sig = xr_stats_sw.copy()
xr_stats_sw_sig = xr_stats_sw_sig.where(xr_stats_sw_sig[:,:,3] < 0.01)
c_mask_sw = np.isfinite(xr_stats_sw_sig[:,:,0]).values
c_mask_sw[mask == False] = False



plt.scatter(t2m['t2m_anom'][2568,633,:], cloud['swnet_anom'][2568,633,:])

plt.scatter(t2m['t2m_anom'][2568,633,:], cloud['swnet_diff'][2568,633,:])












