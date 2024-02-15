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
import matplotlib.pyplot as plt
import glob

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

# Define temperature files
files = sorted(glob.glob(path + 'merra-t2m-resample/*'))

#%%

# Import data
data = xr.open_dataset(path + 'final-forcing-grids.nc')
t2m = xr.open_dataset(path + 'final-temp-grids.nc')


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
xr_stats_snowline = xr.apply_ufunc(new_linregress, t2m['t2m'], 
                          data['snowline'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_snow = xr.apply_ufunc(new_linregress, t2m['t2m'], 
                          data['snow'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_ice = xr.apply_ufunc(new_linregress, t2m['t2m'], 
                          data['ice'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

xr_stats_bulk = xr.apply_ufunc(new_linregress, t2m['t2m'], 
                          data['bulk'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

#%%
# Get only signficant (p<0.05) relationships




#%%
# Bin forcing into elevations
elevations = np.arange(0, 3400, 200)

ice, snowline, snow, bulk = [], [], [], []
area = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    snowline.append(np.nanmean(xr_stats_snowline[:,:,0].values[elevation_mask]))
    ice.append(np.nanmean(xr_stats_ice[:,:,0].values[elevation_mask]))
    snow.append(np.nanmean(xr_stats_snow[:,:,0].values[elevation_mask]))
    bulk.append(np.nanmean(xr_stats_bulk[:,:,0].values[elevation_mask]))

    
#%%

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 6))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(ice, elevations[:-1], color=c1, zorder=2, alpha=0.8, label='Glacier ice')
ax1.plot(snowline, elevations[:-1], color=c2, zorder=2, alpha=0.8, label='Snowline')
ax1.plot(snow, elevations[:-1], color=c3, zorder=2, alpha=0.8, label='Snow')
ax1.plot(bulk, elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Bulk')

ax1.set_ylim(0, 3200)
ax1.legend()

ax2.barh(range(len(area)), area, align='edge',  alpha=0.2, color='blue', edgecolor='k')
ax2.set_ylim(0,17)
ax2.set_yticklabels([])
ax2.set_yticks([])

#%%

# Compute feedback strength at JAR1 69.498, -49.681
#lat, lon = 69.498, -49.681
lat, lon = 67.0955, -49.951 # KAN-L


# First, find the index of the grid point nearest a specific lat/lon.   
abslat = np.abs(data['latitude'] - lat)
abslon = np.abs(data['longitude'] - lon)
c = np.maximum(abslon, abslat)

([xloc], [yloc]) = np.where(c == np.min(c))

print(xr_stats_snowline[xloc,yloc,0])
print(xr_stats_snow[xloc,yloc,0])
print(xr_stats_ice[xloc,yloc,0])

#%%

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'empirical_albedo_model.nc', 
                          'w', format='NETCDF4_CLASSIC')
print('Creating... %s' % path + 'empirical_albedo_model.nc')
dataset.Title = "Slopes and intercepts for temperature vs. albedo relationship"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C., Smith. L. C., Cooley, S. W., and Pearson, B. (in review), Emerging importance of clouds for Greenland Ice Sheet energy balance and meltwater production."
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', linear_inter.shape[0])
lon_dim = dataset.createDimension('x', linear_inter.shape[1])

# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
slope_nc = dataset.createVariable('slope', np.float32, ('y','x'))
inter_nc = dataset.createVariable('intercept', np.float32, ('y','x'))
temp_mean_nc = dataset.createVariable('mean_temp', np.float32, ('y','x'))

# Write data to layers
Y[:] = ds['latitude'].values
X[:] = ds['longitude'].values
slope_nc[:] = linear_slope
inter_nc[:] = linear_inter
temp_mean_nc[:] = np.mean(ds['t2m'], axis=2)

print('Writing data to %s' % path + 'empirical_albedo_model.nc')
    
# Close dataset
dataset.close()




















