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

# Compute anomalies
data['swnet_anom'] = data['swnet'] - data['swnet'].mean(dim='z')
t2m['t2m_anom'] = t2m['t2m'] - t2m['t2m'].mean(dim='z')
swd['swd_allsky_anom'] = swd['swd'] - swd['swd'].mean(dim='z')

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

#%%
xr_stats_cloud = xr.apply_ufunc(new_linregress, swd['swd_allsky_anom'], 
                          data['swnet_anom'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

#%%
# Get only signficant (p<0.05) relationships
xr_stats_temp_sig = xr_stats_temp.copy()
xr_stats_temp_sig = xr_stats_temp_sig.where(xr_stats_temp_sig[:,:,3] < 0.1)

xr_stats_cloud_sig = xr_stats_cloud.copy()
xr_stats_cloud_sig = xr_stats_cloud_sig.where(xr_stats_cloud_sig[:,:,3] < 0.1)

t_mask = np.isfinite(xr_stats_temp_sig[:,:,0])
c_mask = np.isfinite(xr_stats_cloud_sig[:,:,0])

#%%

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

xr_stats_bulk = xr.apply_ufunc(new_linregress, t2m['t2m_anom'], 
                          data['bulk'],
                          input_core_dims=[['z'], ['z']],
                          output_core_dims=[["parameter"]],
                          vectorize=True,
                          dask="parallelized",
                          output_dtypes=['float64'],
                          output_sizes={"parameter": 5})

#%%




#%%

# Bin forcing into elevations
elevations = np.arange(0, 3400, 200)

ice, snowline, snow = [], [], []
area_all, area_valid_temp, area_valid_snow = [], [], []
bulk_t = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area_all.append(elevation_mask.sum())
    area_valid_t = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (t_mask == True)
    area_valid_temp.append(area_valid_t.values.sum())
    
    snowline.append(np.nanmean(xr_stats_snowline[:,:,0].values[area_valid_t]))
    ice.append(np.nanmean(xr_stats_ice[:,:,0].values[area_valid_t]))
    bulk_t.append(np.nanmean(xr_stats_temp_sig[:,:,0].values[area_valid_t]))
    
    area_valid_t = area_valid_t.values
    area_valid_t[snow_mask['snow_mask'].values == 0] = False
    snow.append(np.nanmean(xr_stats_snow[:,:,0].values[area_valid_t]))
    area_valid_snow.append(area_valid_t.sum())


area_all = np.array(area_all)
area_valid_temp = np.array(area_valid_temp)
area_valid_snow = np.array(area_valid_snow)

ice = np.array(ice)
snowline = np.array(snowline)
snow = np.array(snow)

snow[area_valid_snow < 7000] = np.nan
ice[area_valid_temp < 7000] = np.nan
snowline[area_valid_temp < 7000] = np.nan

#%%

coeffs['bulk'] = coeffs['runoff'] * (bulk_t/coeffs['watts'])
coeffs['ice'] = coeffs['runoff'] * (ice/coeffs['watts'])
coeffs['snowline'] = coeffs['runoff'] * (snowline/coeffs['watts'])
coeffs['snow'] = coeffs['runoff'] * (snow/coeffs['watts'])


coeffs['runoff_gt'] = coeffs['runoff']/1e+06*area_all*92
coeffs['runoff_gt'] = coeffs['runoff']/1e+06*area_all*92


coeffs['bulk_gt'] = coeffs['runoff_gt'] * (bulk_t/coeffs['watts'])
coeffs['ice_gt'] = coeffs['runoff_gt'] * (ice/coeffs['watts'])
coeffs['snowline_gt'] = coeffs['runoff_gt'] * (snowline/coeffs['watts'])
coeffs['snow_gt'] = coeffs['runoff_gt'] * (snow/coeffs['watts'])


#%%

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(bulk_t, elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Bulk')

ax1.set_ylim(0, 3200)
ax1.legend(fontsize=13)

ax2.barh(range(len(area_all)), area_all, align='edge',  alpha=0.2, color='blue', edgecolor='k')
ax2.barh(range(len(area_valid_temp)), area_valid_temp, align='edge',  alpha=0.2, color=c1, edgecolor='k')
#ax2.barh(range(len(area_valid_snow)), area_valid_snow, align='edge',  alpha=0.2, color=c2, edgecolor='k')

ax2.set_ylim(0,17)
ax2.set_yticklabels([])

ax3.plot(coeffs['bulk_gt'] , elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Bulk')
ax3.set_ylim(0, 3200)
#ax3.legend(fontsize=13)
ax3.set_yticklabels([])

ax4.plot(ice, elevations[:-1], color=c1, lw=2, zorder=2, alpha=0.8, label='Glacier ice')
ax4.plot(snowline, elevations[:-1], color=c2, lw=2, zorder=2, alpha=0.8, label='Snowline')
ax4.plot(snow, elevations[:-1], color=c3, lw=2, zorder=2, alpha=0.8, label='Snow')
ax4.legend(fontsize=13)
ax4.set_ylim(0, 3200)

ax5.plot(coeffs['ice'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax5.plot(coeffs['snowline'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax5.plot(coeffs['snow'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax5.set_ylim(0, 3200)
#ax5.legend(fontsize=13)
ax5.set_yticklabels([])

ax6.plot(coeffs['ice_gt'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax6.plot(coeffs['snowline_gt'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax6.plot(coeffs['snow_gt'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax6.set_ylim(0, 3200)
#ax6.legend(fontsize=13)
ax6.set_yticklabels([])

ax1.set_xlabel('Radiative feedback (W m$^{-2}$ K$^{-1}$)', fontsize=14)
ax2.set_xlabel('Ice sheet area (km$^2$)', fontsize=14)
ax3.set_xlabel('Meltwater runoff (Gt yr$^{-1}$ K$^{-1}$)', fontsize=14)
ax4.set_xlabel('Radiative feedback (W m$^{-2}$ K$^{-1}$)', fontsize=14)
ax5.set_xlabel('Meltwater runoff (mm w.e. K$^{-1}$)', fontsize=14)
ax6.set_xlabel('Meltwater runoff (Gt yr$^{-1}$ K$^{-1}$)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)
ax4.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=13)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.89, "d", fontsize=24, transform=ax4.transAxes)
ax5.text(0.03, 0.89, "e", fontsize=24, transform=ax5.transAxes)
ax6.text(0.03, 0.89, "f", fontsize=24, transform=ax6.transAxes)

#%%

# Compute feedback strength at JAR1 69.498, -49.681
lat, lon = 69.498, -49.681
#lat, lon = 67.0955, -49.951 # KAN-L


# First, find the index of the grid point nearest a specific lat/lon.   
abslat = np.abs(data['latitude'] - lat)
abslon = np.abs(data['longitude'] - lon)
c = np.maximum(abslon, abslat)

([xloc], [yloc]) = np.where(c == np.min(c))

print(xr_stats_temp_sig[xloc,yloc,0])


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




















