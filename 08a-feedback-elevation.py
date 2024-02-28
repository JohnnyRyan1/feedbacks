#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Radiative feedback analysis.

"""

# Import modules
import xarray as xr
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
data = xr.open_dataset(path + 'feedback-stats.nc')
snow_mask = xr.open_dataset(path + 'snow-mask.nc')

#%%

# Get only signficant (p<0.05) relationships
xr_stats_temp_sig = data['bulk_vs_temp'].copy()
xr_stats_temp_sig = xr_stats_temp_sig.where(xr_stats_temp_sig[:,:,3] < 0.1)

xr_stats_cloud_sig = data['bulk_vs_cloud'].copy()
xr_stats_cloud_sig = xr_stats_cloud_sig.where(xr_stats_cloud_sig[:,:,3] < 0.1)

t_mask = np.isfinite(xr_stats_temp_sig[:,:,0])
c_mask = np.isfinite(xr_stats_cloud_sig[:,:,0])

#%%

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

ice, snowline, snow = [], [], []
area_all, area_valid_temp, area_valid_snow = [], [], []
bulk_t, bulk_c = [], []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area_all.append(elevation_mask.sum())
    area_valid_t = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (t_mask == True)
    area_valid_temp.append(area_valid_t.values.sum())
    
    snowline.append(np.nanmean(data['snowline'].values[:,:,0][area_valid_t]))
    ice.append(np.nanmean(data['ice'].values[:,:,0][area_valid_t]))
    bulk_t.append(np.nanmean(data['bulk_vs_temp'].values[:,:,0][area_valid_t]))
    
    area_valid_c = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (c_mask == True)
    bulk_c.append(np.nanmean(data['bulk_vs_cloud'].values[:,:,0][area_valid_c]))
    
    area_valid_t = area_valid_t.values
    area_valid_t[snow_mask['snow_mask'].values == 0] = False
    snow.append(np.nanmean(data['snow'].values[:,:,0][area_valid_t]))
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

coeffs['bulk'] = coeffs['runoff'] * (bulk_t / coeffs['watts'])
coeffs['ice'] = coeffs['runoff'] * (ice / coeffs['watts'])
coeffs['snowline'] = coeffs['runoff'] * (snowline / coeffs['watts'])
coeffs['snow'] = coeffs['runoff'] * (snow/coeffs['watts'])


coeffs['runoff_gt'] = coeffs['runoff'] / 1e+06 * area_all * 92
coeffs['runoff_gt'] = coeffs['runoff'] / 1e+06 * area_all * 92


coeffs['bulk_gt'] = coeffs['runoff_gt'] * (bulk_t / coeffs['watts'])
coeffs['ice_gt'] = coeffs['runoff_gt'] * (ice / coeffs['watts'])
coeffs['snowline_gt'] = coeffs['runoff_gt'] * (snowline / coeffs['watts'])
coeffs['snow_gt'] = coeffs['runoff_gt'] * (snow / coeffs['watts'])

#%%


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(bulk_t, elevations[:-1], color=c4, zorder=2, alpha=0.8, label='')
ax1.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax1.set_ylim(0, 3400)
ax1.legend(fontsize=13)

ax2.barh(range(len(area_all)), area_all, align='edge',  alpha=0.2, color='blue', edgecolor='k')
ax2.barh(range(len(area_valid_temp)), area_valid_temp, align='edge',  alpha=0.2, color=c1, edgecolor='k')
#ax2.barh(range(len(area_valid_snow)), area_valid_snow, align='edge',  alpha=0.2, color=c2, edgecolor='k')

ax2.set_ylim(0,17)
ax2.set_yticklabels([])

ax3.plot(coeffs['bulk_gt'] , elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Bulk')
ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_ylim(0, 3400)
#ax3.legend(fontsize=13)
ax3.set_yticklabels([])

ax4.plot(ice, elevations[:-1], color=c1, lw=2, zorder=2, alpha=0.8, label='Glacier ice')
ax4.plot(snowline, elevations[:-1], color=c2, lw=2, zorder=2, alpha=0.8, label='Snowline')
ax4.plot(snow, elevations[:-1], color=c3, lw=2, zorder=2, alpha=0.8, label='Snow')
ax4.legend(fontsize=13)
ax4.set_ylim(0, 3400)

ax5.plot(coeffs['ice'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax5.plot(coeffs['snowline'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax5.plot(coeffs['snow'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax5.set_ylim(0, 3400)
#ax5.legend(fontsize=13)
ax5.set_yticklabels([])

ax6.plot(coeffs['ice_gt'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax6.plot(coeffs['snowline_gt'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax6.plot(coeffs['snow_gt'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax6.set_ylim(0, 3400)
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

# Fraction of significant correlations between air temp and SWnet <1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values < 1600)
print(t_mask.values[elevation_mask].sum() / mask[elevation_mask].sum())

# Mean slope of SWnet vs. air temp <1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values < 1600) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))

# Mean slope of SWnet vs. air temp 0-200 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values < 200) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))

# Mean slope of SWnet vs. air temp 1400-1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 1400) & (ismip_1km['SRF'].values < 1600) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))

# Percentage increase in radiative forcing 
sig_mask = (mask == True) & (ismip_1km['SRF'].values < 1600) & (t_mask == True)
print((np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask])+38.1)/38.1)

# Fraction of significant correlations between air temp and SWnet >1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 1600) & (t_mask == True)
ele_mask = (mask == True) & (ismip_1km['SRF'].values > 1600)
print(sig_mask.sum().values/ele_mask.sum())

# Fraction significant correlations between air temp and SWnet for entire ice sheet
print(t_mask.sum().values/mask.sum())















