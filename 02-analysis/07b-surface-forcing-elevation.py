#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Surface radiative forcing analysis.

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Library/CloudStorage/OneDrive-DukeUniversity/research/feedbacks/re-revision/'

#%%

# Import data
data = xr.open_dataset(path + 'final-surface-forcing-grids.nc')

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-melt-coeffs.csv')
factor = pd.read_csv(path + 'melt-factors.csv')

#%%

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

ice, snowline, snow, swnet = [], [], [], []
ice_std, snowline_std, snow_std, swnet_std = [], [], [], []
area = []

for e in range(len(elevations) - 1):
    print(e)
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    ice.append(np.nanmean(np.nanmean(np.nanmean(data['ice'].values, axis=2)[elevation_mask])))
    snowline.append(np.nanmean(np.nanmean(data['snowline'].values, axis=2)[elevation_mask]))
    snow.append(np.nanmean(np.nanmean(data['snow'].values, axis=2)[elevation_mask]))
    swnet.append(np.nanmean(np.nanmean(data['swnet'].values, axis=2)[elevation_mask]))
    
    ice_std.append(np.nanmean(np.nanmean(np.nanstd(data['ice'].values, axis=2)[elevation_mask])))
    snowline_std.append(np.nanmean(np.nanstd(data['snowline'].values, axis=2)[elevation_mask]))
    snow_std.append(np.nanmean(np.nanstd(data['snow'].values, axis=2)[elevation_mask]))
    swnet_std.append(np.nanmean(np.nanstd(data['swnet'].values, axis=2)[elevation_mask]))

area = np.array(area)
ice = np.array(ice)
snowline = np.array(snowline)
snow = np.array(snow)

ice_std = np.array(ice_std)
snowline_std = np.array(snowline_std)
snow_std = np.array(snow_std)

# Save as DataFrame
df = pd.DataFrame(list(zip(ice, snowline, snow)))
df.columns = ['ice', 'snowline', 'snow']
df['bulk'] = df['ice'] + df['snowline'] + df['snow']
#df.to_csv(path + 'radiative-forcing-elevation.csv', index=False)

#%%

# Compute melt factor
coeffs['factor'] = coeffs['melt'] / coeffs['watts']

coeffs['ice'] = coeffs['factor'] * ice
coeffs['snowline'] = coeffs['factor'] * snowline
coeffs['snow'] = coeffs['factor'] * snow

coeffs['ice_std'] = coeffs['factor'] * ice_std
coeffs['snowline_std'] = coeffs['factor'] * snowline_std
coeffs['snow_std'] = coeffs['factor'] * snow_std

# Convert to Gt (mm to km, multiply by area, multiply by number of days)
coeffs['melt_gt'] = coeffs['melt']/1e+06*area*92

coeffs['ice_gt'] = coeffs['ice']/1e+06*area*92
coeffs['snowline_gt'] = coeffs['snowline']/1e+06*area*92
coeffs['snow_gt'] = coeffs['snow']/1e+06*area*92

coeffs['ice_gt_std'] = coeffs['ice_std']/1e+06*area*92
coeffs['snowline_gt_std'] = coeffs['snowline_std']/1e+06*area*92
coeffs['snow_gt_std'] = coeffs['snow_std']/1e+06*area*92

#coeffs.to_csv(path + 'final-coeffs.csv')
#%%

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(ice+snowline+snow, elevations[:-1], color=c4, zorder=2, lw=2, 
         alpha=0.8, label='Surface radiative forcing')
ax1.fill_betweenx(elevations[:-1],
                 ice+snowline+snow+ice_std+snowline_std+snow_std,
                 ice+snowline+snow-(ice_std+snowline_std+snow_std),
                 zorder=1, color=c4, alpha=0.2)
ax1.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax1.legend(fontsize=13)
ax1.set_xlim(0, 100)

ax2.barh(range(len(area)), area, align='edge',  alpha=0.4, color=c4, edgecolor='k')
ax2.set_ylim(0,17)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=13)
ax2.set_yticklabels([])

#ax3.plot(coeffs['melt_gt'], elevations[:-1], color=c4, zorder=2, 
#         lw=2, alpha=0.8, label='Total')
ax3.plot(coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt'], 
         elevations[:-1], color=c4, zorder=2, 
         lw=2, alpha=0.8, label='')
ax3.fill_betweenx(elevations[:-1],
                 coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt']+
                 coeffs['ice_gt_std']+coeffs['snowline_gt_std']+coeffs['snow_gt_std'],
                 coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt']-
                 (coeffs['ice_gt_std']+coeffs['snowline_gt_std']+coeffs['snow_gt_std']),
                 zorder=1, color=c4, alpha=0.2)
ax3.set_yticklabels([])
ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_xlim(0, 8.5)

ax4.plot(ice, elevations[:-1], color=c1, zorder=2, lw=2, alpha=0.8, label='Glacier ice')
ax4.fill_betweenx(elevations[:-1], ice+ice_std, ice-ice_std, color=c1, zorder=1, alpha=0.2)
ax4.plot(snowline, elevations[:-1], color=c2, lw=2, zorder=2, alpha=0.8, label='Snowline')
ax4.fill_betweenx(elevations[:-1], snowline+snowline_std, snowline-snowline_std, color=c2, zorder=1, alpha=0.2)
ax4.plot(snow, elevations[:-1], color=c3, lw=2, zorder=2, alpha=0.8, label='Snow')
ax4.fill_betweenx(elevations[:-1], snow+snow_std, snow-snow_std, color=c3, zorder=1, alpha=0.2)
ax4.set_xlim(0, 60)
ax4.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax4.legend(fontsize=13)

ax5.plot(coeffs['ice'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax5.fill_betweenx(elevations[:-1], coeffs['ice']+coeffs['ice_std'], 
                  coeffs['ice']-coeffs['ice_std'], color=c1, zorder=1, alpha=0.2)
ax5.plot(coeffs['snowline'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax5.fill_betweenx(elevations[:-1], coeffs['snowline']+coeffs['snowline_std'], 
                  coeffs['snowline']-coeffs['snowline_std'], color=c2, zorder=1, alpha=0.2)
ax5.plot(coeffs['snow'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax5.fill_betweenx(elevations[:-1], coeffs['snow']+coeffs['snow_std'], 
                  coeffs['snow']-coeffs['snow_std'], color=c3, zorder=1, alpha=0.2)
ax5.set_xlim(0, 2.3)
ax5.set_yticklabels([])
ax5.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)

ax6.plot(coeffs['ice_gt'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax6.fill_betweenx(elevations[:-1], coeffs['ice_gt']+coeffs['ice_gt_std'], 
                  coeffs['ice_gt']-coeffs['ice_gt_std'], color=c1, zorder=1, alpha=0.2)
ax6.plot(coeffs['snowline_gt'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax6.fill_betweenx(elevations[:-1], coeffs['snowline_gt']+coeffs['snowline_gt_std'], 
                  coeffs['snowline_gt']-coeffs['snowline_gt_std'], color=c2, zorder=1, alpha=0.2)
ax6.plot(coeffs['snow_gt'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax6.fill_betweenx(elevations[:-1], coeffs['snow_gt']+coeffs['snow_gt_std'], 
                  coeffs['snow_gt']-coeffs['snow_gt_std'], color=c3, zorder=1, alpha=0.2)
ax6.set_xlim(0, 4)
ax6.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax6.set_yticklabels([])


ax1.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax2.set_xlabel('Ice sheet area (km$^2$)', fontsize=14)
ax3.set_xlabel('Meltwater production (Gt yr$^{-1}$)', fontsize=14)
ax4.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax5.set_xlabel('Meltwater production (mm w.e. d$^{-1}$)', fontsize=14)
ax6.set_xlabel('Meltwater production (Gt yr$^{-1}$)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)
ax4.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

for ax in [ax1, ax3, ax4, ax5, ax6]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_ylim(0, 3400)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.89, "d", fontsize=24, transform=ax4.transAxes)
ax5.text(0.03, 0.89, "e", fontsize=24, transform=ax5.transAxes)
ax6.text(0.03, 0.89, "f", fontsize=24, transform=ax6.transAxes)

ax1.set_xlim(0, 100)
ax4.set_xlim(0, 100)
ax3.set_xlim(0, 8.6)
ax6.set_xlim(0, 8.6)

fig.savefig(savepath + 'fig2.png', dpi=300)


#%%

"""
P2
"""

# Mean radiative forcing between 0-1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 1600)
print(np.mean((np.mean(data['ice'], axis=2) + \
               np.mean(data['snowline'], axis=2) + \
               np.mean(data['snow'], axis=2)).values[elevation_mask]))
print(np.mean((np.std(data['ice'], axis=2) + \
               np.std(data['snowline'], axis=2) + \
               np.std(data['snow'], axis=2)).values[elevation_mask]))

    
# Mean radiative forcing between >1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 1600)
print(np.mean((np.mean(data['ice'], axis=2) + \
               np.mean(data['snowline'], axis=2) + \
               np.mean(data['snow'], axis=2)).values[elevation_mask]))
print(np.mean((np.std(data['ice'], axis=2) + \
               np.std(data['snowline'], axis=2) + \
               np.std(data['snow'], axis=2)).values[elevation_mask]))
    
#%%
"""
P3
"""
# Fraction of melt produced between 0-1600 m
print(np.sum(coeffs['melt_gt'][0:8]) / np.sum(coeffs['melt_gt']))

# Mean radiative forcing from snowlines between 0-1600 m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 1600)
print(np.mean(np.mean(data['snowline'], axis=2).values[elevation_mask]))
print(np.mean(np.std(data['snowline'], axis=2).values[elevation_mask]))

# Mean radiative forcing from glacier ice between 0-1600 m
print(np.mean(np.mean(data['ice'], axis=2).values[elevation_mask]))
print(np.mean(np.std(data['ice'], axis=2).values[elevation_mask]))

# Mean radiative forcing from snow between 0-1600 m
print(np.mean(np.mean(data['snow'], axis=2).values[elevation_mask]))
print(np.mean(np.std(data['snow'], axis=2).values[elevation_mask]))

#%%
"""
P5
"""

# Amount of meltwater melt due to snowline RF
print(np.sum(coeffs['snowline_gt']), '+/-', np.sum(coeffs['snowline_gt_std']))

# Amount of meltwater melt due to snowline RF compared with total melt
print(np.sum(coeffs['snowline_gt']) / np.sum(coeffs['melt_gt']))

# Amount of meltwater melt due to snowline RF compared with total melt
print((np.sum(coeffs['snowline_gt']) + np.sum(coeffs['snowline_gt_std'])) / np.sum(coeffs['melt_gt']))

# Amount of meltwater melt due to snow RF
print(np.sum(coeffs['snow_gt']), '+/-', np.sum(coeffs['snow_gt_std']))

# Amount of meltwater melt due to ice RF compared with total melt
print(np.sum(coeffs['snow_gt']) / np.sum(coeffs['melt_gt']))

# Amount of meltwater melt due to ice RF compared with total melt
print((np.sum(coeffs['snow_gt']) + np.sum(coeffs['snow_gt_std'])) / np.sum(coeffs['melt_gt']))

# Amount of meltwater melt due to ice RF
print(np.sum(coeffs['ice_gt']), '+/-', np.sum(coeffs['ice_gt_std']))

# Amount of meltwater melt due to ice RF compared with total melt
print(np.sum(coeffs['ice_gt']) / np.sum(coeffs['melt_gt']))

# Amount of meltwater melt due to ice RF compared with total melt
print((np.sum(coeffs['ice_gt']) + np.sum(coeffs['ice_gt_std'])) / np.sum(coeffs['melt_gt']))


all_gt = np.sum(coeffs['ice_gt']) + np.sum(coeffs['snowline_gt']) + np.sum(coeffs['snow_gt'])

all_gt_max = np.sum(coeffs['ice_gt']) + np.sum(coeffs['snowline_gt']) + np.sum(coeffs['snow_gt']) +\
         np.sum(coeffs['ice_gt_std']) + np.sum(coeffs['snowline_gt_std']) + np.sum(coeffs['snow_gt_std'])


# Area difference between 0-800 m and 800-1600 m
print((area[4:8].sum()-area[0:4].sum())/area[0:4].sum())

# Contribution of RF to total melt
print(all_gt / np.sum(coeffs['melt_gt']))

# Contribution of RF to total melt
print(all_gt_max / np.sum(coeffs['melt_gt']))


# Amount of meltwater melt due to surface RF
print(all_gt, '+/-', all_gt_max - all_gt)
               

#%%

"""
Uncertainty due to melt factors
"""

# Compute melt factor
coeffs['factor'] = coeffs['melt'] / coeffs['watts']
coeffs['factor_max'] = coeffs['factor'] + factor.std(axis=1)

# Mean
coeffs['ice'] = coeffs['factor'] * ice
coeffs['snowline'] = coeffs['factor'] * snowline
coeffs['snow'] = coeffs['factor'] * snow

# Convert to Gt (mm to km, multiply by area, multiply by number of days)
coeffs['melt_gt'] = coeffs['melt']/1e+06*area*92

coeffs['ice_gt'] = coeffs['ice']/1e+06*area*92
coeffs['snowline_gt'] = coeffs['snowline']/1e+06*area*92
coeffs['snow_gt'] = coeffs['snow']/1e+06*area*92

# Max
coeffs['ice_max'] = coeffs['factor_max'] * ice
coeffs['snowline_max'] = coeffs['factor_max'] * snowline
coeffs['snow_max'] = coeffs['factor_max'] * snow

coeffs['ice_gt_max'] = coeffs['ice_max']/1e+06*area*92
coeffs['snowline_gt_max'] = coeffs['snowline_max']/1e+06*area*92
coeffs['snow_gt_max'] = coeffs['snow_max']/1e+06*area*92

               
all_gt = np.sum(coeffs['ice_gt']) + np.sum(coeffs['snowline_gt']) + np.sum(coeffs['snow_gt'])
all_gt_max = np.sum(coeffs['ice_gt_max']) + np.sum(coeffs['snowline_gt_max']) + np.sum(coeffs['snow_gt_max'])


#%%

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(ice+snowline+snow, elevations[:-1], color=c4, zorder=2, lw=2, 
         alpha=0.8, label='Surface radiative forcing')
ax1.fill_betweenx(elevations[:-1],
                 ice+snowline+snow+ice_std+snowline_std+snow_std,
                 ice+snowline+snow-(ice_std+snowline_std+snow_std),
                 zorder=1, color=c4, alpha=0.2)
ax1.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax1.legend(fontsize=13)
ax1.set_xlim(0, 100)

ax2.barh(range(len(area)), area, align='edge',  alpha=0.4, color=c4, edgecolor='k')
ax2.set_ylim(0,17)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=13)
ax2.set_yticklabels([])

#ax3.plot(coeffs['melt_gt'], elevations[:-1], color=c4, zorder=2, 
#         lw=2, alpha=0.8, label='Total')
ax3.plot(coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt'], 
         elevations[:-1], color=c4, zorder=2, 
         lw=2, alpha=0.8, label='')
ax3.fill_betweenx(elevations[:-1],
                 coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt']+
                 coeffs['ice_gt_std']+coeffs['snowline_gt_std']+coeffs['snow_gt_std'],
                 coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt']-
                 (coeffs['ice_gt_std']+coeffs['snowline_gt_std']+coeffs['snow_gt_std']),
                 zorder=1, color=c4, alpha=0.2)
ax3.legend(fontsize=13)
ax3.set_yticklabels([])
ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_xlim(0, 8.5)

ax4.plot(snowline, elevations[:-1], color=c2, lw=2, zorder=2, alpha=0.8, label='Snowline')
ax4.fill_betweenx(elevations[:-1], snowline+snowline_std, snowline-snowline_std, color=c2, zorder=1, alpha=0.2)
ax4.set_xlim(0, 60)
ax4.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax4.legend(fontsize=13)

ax5.plot(coeffs['snowline'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax5.fill_betweenx(elevations[:-1], coeffs['snowline']+coeffs['snowline_std'], 
                  coeffs['snowline']-coeffs['snowline_std'], color=c2, zorder=1, alpha=0.2)
ax5.set_xlim(0, 2.3)
ax5.set_yticklabels([])
ax5.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)

ax6.plot(coeffs['snowline_gt'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax6.fill_betweenx(elevations[:-1], coeffs['snowline_gt']+coeffs['snowline_gt_std'], 
                  coeffs['snowline_gt']-coeffs['snowline_gt_std'], color=c2, zorder=1, alpha=0.2)
ax6.set_xlim(0, 4)
ax6.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax6.set_yticklabels([])


ax1.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax2.set_xlabel('Ice sheet area (km$^2$)', fontsize=14)
ax3.set_xlabel('Meltwater melt (Gt yr$^{-1}$)', fontsize=14)
ax4.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax5.set_xlabel('Meltwater production (mm w.e. d$^{-1}$)', fontsize=14)
ax6.set_xlabel('Meltwater production (Gt yr$^{-1}$)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)
ax4.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

for ax in [ax1, ax3, ax4, ax5, ax6]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_ylim(0, 3400)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.89, "d", fontsize=24, transform=ax4.transAxes)
ax5.text(0.03, 0.89, "e", fontsize=24, transform=ax5.transAxes)
ax6.text(0.03, 0.89, "f", fontsize=24, transform=ax6.transAxes)

#fig.savefig(savepath + 'radiative-forcing-snowline.png', dpi=200)





