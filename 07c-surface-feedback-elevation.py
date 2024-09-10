#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Surface radiative feedback analysis.

"""

# Import modules
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/feedbacks/re-revision/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-melt-coeffs.csv')

# Read radiative forcing
comp = pd.read_csv(path + 'radiative-forcing-elevation.csv')

#%%

# Import data
data = xr.open_dataset(path + 'feedback-stats.nc')
snow_mask = xr.open_dataset(path + 'snow-mask.nc')

#%%

# Get only signficant (p<0.1) relationships
xr_stats_temp_sig = data['surface_feedback'].copy()
xr_stats_temp_sig = xr_stats_temp_sig.where(xr_stats_temp_sig[:,:,3] < 0.1)
t_mask = np.isfinite(xr_stats_temp_sig[:,:,0]).values
t_mask[mask == False] = False

std_err = data['surface_feedback'][:,:,4]

#%%

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

ice, snowline, snow = [], [], []
area_all, area_valid_temp, area_valid_snow = [], [], []
bulk_t, bulk_t_std = [], []
snowline_std, ice_std, snow_std = [], [], []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area_all.append(elevation_mask.sum())
    area_valid_t = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (t_mask == True)
    area_valid_temp.append(area_valid_t.sum())
    
    snowline.append(np.nanmean(data['snowline_feedback'].values[:,:,0][area_valid_t]))
    snowline_std.append(np.nanmean(data['snowline_feedback'].values[:,:,4][area_valid_t]))

    ice.append(np.nanmean(data['ice_feedback'].values[:,:,0][area_valid_t]))
    ice_std.append(np.nanmean(data['ice_feedback'].values[:,:,4][area_valid_t]))

    
    bulk_t.append(np.nanmean(data['surface_feedback'].values[:,:,0][area_valid_t]))
    bulk_t_std.append(np.nanmean(data['surface_feedback'].values[:,:,4][area_valid_t]))

    area_valid_t[snow_mask['snow_mask'].values == 0] = False
    snow.append(np.nanmean(data['snow_feedback'].values[:,:,0][area_valid_t]))
    snow_std.append(np.nanmean(data['snow_feedback'].values[:,:,4][area_valid_t]))
    area_valid_snow.append(area_valid_t.sum())
    

area_all = np.array(area_all)
area_valid_temp = np.array(area_valid_temp)
area_valid_snow = np.array(area_valid_snow)

ice = np.array(ice)
snowline = np.array(snowline)
snow = np.array(snow)

ice_std = np.array(ice_std)
snowline_std = np.array(snowline_std)
snow_std = np.array(snow_std)

snow[area_valid_snow < 7000] = np.nan
ice[area_valid_temp < 7000] = np.nan
snowline[area_valid_temp < 7000] = np.nan

#%%
# Compute individual components of feedback
bulk_t_ice = bulk_t * (comp['ice'] / comp['bulk'])
bulk_t_snowline = bulk_t * (comp['snowline'] / comp['bulk'])
bulk_t_snow = bulk_t * (comp['snow'] / comp['bulk'])

bulk_t_ice_std = bulk_t_std * (comp['ice'] / comp['bulk'])
bulk_t_snowline_std = bulk_t_std * (comp['snowline'] / comp['bulk'])
bulk_t_snow_std = bulk_t_std * (comp['snow'] / comp['bulk'])

bulk_t = np.array(bulk_t)
bulk_t_std = np.array(bulk_t_std)

#%%

# Compute melt factor
coeffs['factor'] = coeffs['melt'] / coeffs['watts']

coeffs['ice'] = coeffs['factor'] * bulk_t_ice
coeffs['ice_std'] = coeffs['factor'] * bulk_t_ice_std
coeffs['snowline'] = coeffs['factor'] * bulk_t_snowline
coeffs['snowline_std'] = coeffs['factor'] * bulk_t_snowline_std
coeffs['snow'] = coeffs['factor'] * bulk_t_snow
coeffs['snow_std'] = coeffs['factor'] * bulk_t_snow_std
coeffs['bulk'] = coeffs['factor'] * bulk_t
coeffs['bulk_std'] = coeffs['factor'] * bulk_t_std

coeffs['melt_gt'] = coeffs['melt'] / 1e+06 * area_all * 92

coeffs['ice_gt'] = coeffs['ice']/1e+06*area_all*92
coeffs['ice_gt_std'] = coeffs['ice_std']/1e+06*area_all*92
coeffs['snowline_gt'] = coeffs['snowline']/1e+06*area_all*92
coeffs['snowline_gt_std'] = coeffs['snowline_std']/1e+06*area_all*92
coeffs['snow_gt'] = coeffs['snow']/1e+06*area_all*92
coeffs['snow_gt_std'] = coeffs['snow_std']/1e+06*area_all*92
coeffs['bulk_gt'] = coeffs['bulk']/1e+06*area_all*92
coeffs['bulk_gt_std'] = coeffs['bulk_std']/1e+06*area_all*92


#%%


fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(bulk_t, elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Surface radiative feedback')
ax1.fill_betweenx(elevations[:-1], bulk_t+(bulk_t_std), bulk_t-(bulk_t_std),
                 zorder=1, color=c4, alpha=0.2)
ax1.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax1.set_ylim(0, 3400)
ax1.legend(fontsize=13)
ax1.set_xlim(0,11)

ax2.barh(range(len(area_all)), area_all, align='edge',  alpha=0.2, color='blue', 
         edgecolor='k', label='all')
ax2.barh(range(len(area_valid_temp)), area_valid_temp, align='edge',  alpha=0.2, 
         color=c1, edgecolor='k', label='p<0.1')
ax2.legend(fontsize=13)
ax2.set_ylim(0,17)
ax2.set_yticklabels([])

ax3.plot(coeffs['bulk_gt'], elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Bulk')
ax3.fill_betweenx(elevations[:-1], coeffs['bulk_gt']+coeffs['bulk_gt_std'], 
                  coeffs['bulk_gt']-coeffs['bulk_gt_std'],
                  zorder=1, color=c4, alpha=0.2)

ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_ylim(0, 3400)
#ax3.legend(fontsize=13)
ax3.set_yticklabels([])

ax3.plot(coeffs['bulk_gt'] , elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Bulk')
ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_ylim(0, 3400)
ax3.set_xlim(0,1.6)
ax3.set_yticklabels([])

ax4.plot(bulk_t_ice, elevations[:-1], color=c1, lw=2, zorder=2, alpha=0.8, label='Glacier ice')
ax4.fill_betweenx(elevations[:-1], bulk_t_ice+bulk_t_ice_std, bulk_t_ice-bulk_t_ice_std,
                 zorder=1, color=c1, alpha=0.2)
ax4.plot(bulk_t_snowline, elevations[:-1], color=c2, lw=2, zorder=2, alpha=0.8, label='Snowline')
ax4.fill_betweenx(elevations[:-1], bulk_t_snowline+bulk_t_snowline_std, 
                  bulk_t_snowline-bulk_t_snowline_std,
                 zorder=1, color=c2, alpha=0.2)
ax4.plot(bulk_t_snow, elevations[:-1], color=c3, lw=2, zorder=2, alpha=0.8, label='Snow')
ax4.fill_betweenx(elevations[:-1], bulk_t_snow+bulk_t_snow_std, 
                  bulk_t_snow-bulk_t_snow_std,
                 zorder=1, color=c3, alpha=0.2)
ax4.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax4.legend(fontsize=13)
ax4.set_ylim(0, 3400)
ax4.set_xlim(0,8)

ax5.plot(coeffs['ice'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax5.fill_betweenx(elevations[:-1], coeffs['ice']+coeffs['ice_std'], 
                  coeffs['ice']-coeffs['ice_std'],
                  zorder=1, color=c1, alpha=0.2)
ax5.plot(coeffs['snowline'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax5.fill_betweenx(elevations[:-1], coeffs['snowline']+coeffs['snowline_std'], 
                  coeffs['snowline']-coeffs['snowline_std'],
                  zorder=1, color=c2, alpha=0.2)
ax5.plot(coeffs['snow'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax5.fill_betweenx(elevations[:-1], coeffs['snow']+coeffs['snow_std'], 
                  coeffs['snow']-coeffs['snow_std'],
                  zorder=1, color=c3, alpha=0.2)
ax5.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax5.set_ylim(0, 3400)
ax5.set_xlim(0,0.3)
ax5.set_yticklabels([])

ax6.plot(coeffs['ice_gt'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax6.fill_betweenx(elevations[:-1], coeffs['ice_gt']+coeffs['ice_gt_std'], 
                  coeffs['ice_gt']-coeffs['ice_gt_std'],
                  zorder=1, color=c1, alpha=0.2)
ax6.plot(coeffs['snowline_gt'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax6.fill_betweenx(elevations[:-1], coeffs['snowline_gt']+coeffs['snowline_gt_std'], 
                  coeffs['snowline_gt']-coeffs['snowline_gt_std'],
                  zorder=1, color=c2, alpha=0.2)
ax6.plot(coeffs['snow_gt'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax6.fill_betweenx(elevations[:-1], coeffs['snow_gt']+coeffs['snow_gt_std'], 
                  coeffs['snow_gt']-coeffs['snow_gt_std'],
                  zorder=1, color=c3, alpha=0.2)
ax6.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax6.set_ylim(0, 3400)
ax6.set_xlim(0,1)
ax6.set_yticklabels([])

ax1.set_xlabel('Radiative feedback (W m$^{-2}$ K$^{-1}$)', fontsize=14)
ax2.set_xlabel('Ice sheet area (km$^2$)', fontsize=14)
ax3.set_xlabel('Meltwater production (Gt yr$^{-1}$ K$^{-1}$)', fontsize=14)
ax4.set_xlabel('Radiative feedback (W m$^{-2}$ K$^{-1}$)', fontsize=14)
ax5.set_xlabel('Meltwater production (mm w.e. K$^{-1}$)', fontsize=14)
ax6.set_xlabel('Meltwater production (Gt yr$^{-1}$ K$^{-1}$)', fontsize=14)

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

ax1.set_xlim(0, 11.5)
ax4.set_xlim(0, 11.5)
ax3.set_xlim(0, 1.6)
ax6.set_xlim(0, 1.6)

fig.savefig(savepath + 'fig3.png', dpi=300)


#%%

# Fraction of significant correlations between air temp and SWnet <1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values < 1600)
print(t_mask[elevation_mask].sum() / mask[elevation_mask].sum())

# Mean slope of SWnet vs. air temp <1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values < 1600) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))
print(np.nanmean(xr_stats_temp_sig[:,:,4].values[sig_mask]))

# Mean slope of SWnet vs. air temp 0-200 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values < 200) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))
print(np.nanmean(xr_stats_temp_sig[:,:,4].values[sig_mask]))

# Mean slope of SWnet vs. air temp 1400-1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 1400) & (ismip_1km['SRF'].values < 1600) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))
print(np.nanmean(xr_stats_temp_sig[:,:,4].values[sig_mask]))

# Percentage increase in radiative forcing 
sig_mask = (mask == True) & (ismip_1km['SRF'].values < 1600) & (t_mask == True)
print((np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask])+36.8)/36.8)

# Mean slope of SWnet vs. air temp >1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 1600) & (t_mask == True)
print(np.nanmean(xr_stats_temp_sig[:,:,0].values[sig_mask]))
print(np.nanmean(xr_stats_temp_sig[:,:,4].values[sig_mask]))

# Fraction of significant correlations between air temp and SWnet >1600 m
sig_mask = (mask == True) & (ismip_1km['SRF'].values > 1600) & (t_mask == True)
ele_mask = (mask == True) & (ismip_1km['SRF'].values > 1600)
print(sig_mask.sum()/ele_mask.sum())

#%%


print(coeffs['bulk_gt'].sum(), '+/-', coeffs['bulk_gt_std'].sum())
















