#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analysis at the grid cell scale

"""

# Import modules
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

#%%

# Import data
data = xr.open_dataset(path + 'final-forcing-grids.nc')

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-runoff-coeffs.csv')
factor = pd.read_csv(path + 'runoff-factors.csv')

#%%

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

ice, snowline, snow, swnet = [], [], [], []
area = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    ice.append(np.nanmean(data['ice'].values[elevation_mask]))
    snowline.append(np.nanmean(data['snowline'].values[elevation_mask]))
    snow.append(np.nanmean(data['snow'].values[elevation_mask]))
    swnet.append(np.nanmean(data['swnet'].values[elevation_mask]))

area = np.array(area)
ice = np.array(ice)
snowline = np.array(snowline)
snow = np.array(snow)

# Compute runoff factor
coeffs['factor'] = coeffs['runoff'] / coeffs['watts']

coeffs['ice'] = coeffs['factor'] * ice
coeffs['snowline'] = coeffs['factor'] * snowline
coeffs['snow'] = coeffs['factor'] * snow

# Convert to Gt (mm to km, multiply by area, multiply by number of days)
coeffs['runoff_gt'] = coeffs['runoff']/1e+06*area*92

coeffs['ice_gt'] = coeffs['ice']/1e+06*area*92
coeffs['snowline_gt'] = coeffs['snowline']/1e+06*area*92
coeffs['snow_gt'] = coeffs['snow']/1e+06*area*92

#%%

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(16, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

#ax1.plot(swnet, elevations[:-1], color=c4, zorder=2, lw=2, alpha=0.8, 
#         label='SW$_{net}$')
ax1.plot(ice+snowline+snow, elevations[:-1], color=c4, zorder=2, lw=2, 
         alpha=0.8, label='SW$_{net}$ due to RF')
ax1.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax1.legend(fontsize=13)
ax1.set_xlim(0, 100)

ax2.barh(range(len(area)), area, align='edge',  alpha=0.4, color=c4, edgecolor='k')
ax2.set_ylim(0,17)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=13)
ax2.set_yticklabels([])

ax3.plot(coeffs['runoff_gt'], elevations[:-1], color=c4, zorder=2, 
         lw=2, alpha=0.8, label='Total')
#ax1.fill_between()
ax3.legend(fontsize=13)
ax3.set_yticklabels([])
ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_xlim(0, 45)

ax4.plot(ice, elevations[:-1], color=c1, zorder=2, lw=2, alpha=0.8, label='Glacier ice')
ax4.plot(snowline, elevations[:-1], color=c2, lw=2, zorder=2, alpha=0.8, label='Snowline')
ax4.plot(snow, elevations[:-1], color=c3, lw=2, zorder=2, alpha=0.8, label='Snow')
ax4.set_xlim(0, 60)
ax4.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax4.legend(fontsize=13)

ax5.plot(coeffs['ice'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax5.plot(coeffs['snowline'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax5.plot(coeffs['snow'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax5.set_xlim(0, 2.3)
ax5.set_yticklabels([])
ax5.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)

ax6.plot(coeffs['ice_gt'], elevations[:-1], color=c1, zorder=2, 
         lw=2, alpha=0.8, label='Glacier ice')
ax6.plot(coeffs['snowline_gt'], elevations[:-1], color=c2, 
         lw=2, zorder=2, alpha=0.8, label='Snowline')
ax6.plot(coeffs['snow_gt'], elevations[:-1], color=c3, 
         lw=2, zorder=2, alpha=0.8, label='Snow')
ax6.set_xlim(0, 3.2)
ax6.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax6.set_yticklabels([])


ax1.set_xlabel('SW radiative forcing (W m$^{-2}$)', fontsize=14)
ax2.set_xlabel('Ice sheet area (km$^2$)', fontsize=14)
ax3.set_xlabel('Meltwater runoff (Gt yr$^{-1}$)', fontsize=14)
ax4.set_xlabel('SW radiative forcing (W m$^{-2}$)', fontsize=14)
ax5.set_xlabel('Meltwater runoff (mm w.e. d$^{-1}$)', fontsize=14)
ax6.set_xlabel('Meltwater runoff (Gt yr$^{-1}$)', fontsize=14)

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

fig.savefig(savepath + 'radiative-forcing-elevation.png', dpi=200)
#%%

fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(factor.mean(axis=1), elevations[:-1], color=c1, zorder=2, s=100, alpha=0.8, label='')
ax1.plot(factor.mean(axis=1), elevations[:-1], color=c1, zorder=2, lw=2, alpha=0.5, 
         ls='dashed', label='')
ax1.fill_betweenx(elevations[:-1], factor.mean(axis=1) - factor.std(axis=1), 
                  factor.mean(axis=1) + factor.std(axis=1), zorder=1,
                  color='grey', alpha=0.3)
ax1.set_ylim(0, 3200)

ax1.set_xlabel('Runoff factor (mm W$^{-1}$)', fontsize=14)
ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)

fig.savefig(savepath + 'watts-to-runoff.png', dpi=200)

#%%

# Stats

# Mean radiative forcing between 0-1000 m
print(np.mean((ice+snowline+snow)[0:8]))

# Mean radiative forcing between >1500 m
print(np.mean((ice+snowline+snow)[8:]))

# Fraction of runoff produced between 0-1600 m
print(np.sum(coeffs['runoff_gt'][0:8]) / np.sum(coeffs['runoff_gt']))

# Mean radiative forcing from snowlines between 0-1600 m
print(np.mean(snowline[0:8]))

# Mean radiative forcing from glacier ice between 0-1600 m
print(np.mean(ice[0:8]))

# Mean radiative forcing from snow between 0-1600 m
print(np.mean(snow[0:8]))

# Contribution of snowline to RF between 0-1600 m
print(np.sum(snowline[0:8]) / np.sum((ice+snowline+snow)[0:8]))

# Contribution of snow to RF between 0-1600 m
print(np.sum(snow[0:8]) / np.sum((ice+snowline+snow)[0:8]))

# Contribution of glacier ice to RF between 0-1600 m
print(np.sum(ice[0:8]) / np.sum((ice+snowline+snow)[0:8]))

# Amount of meltwater runoff due to ice RF
print(np.sum(coeffs['ice_gt']))

# Amount of meltwater runoff due to snowline RF
print(np.sum(coeffs['snowline_gt']))

# Amount of meltwater runoff due to snow RF
print(np.sum(coeffs['snow_gt']))

all_gt = np.sum(coeffs['ice_gt']) + np.sum(coeffs['snowline_gt']) + np.sum(coeffs['snow_gt'])

# Contribution of meltwater runoff due to ice RF compared with total runoff
print(np.sum(coeffs['ice_gt']) / np.sum(coeffs['runoff_gt']))

# Amount of meltwater runoff due to snowline RF compared with total runoff
print(np.sum(coeffs['snowline_gt']) / np.sum(coeffs['runoff_gt']))

# Amount of meltwater runoff due to snow RF compared with total runoff
print(np.sum(coeffs['snow_gt']) / np.sum(coeffs['runoff_gt']))

# Area difference between 0-800 m and 800-1600 m
print((area[4:8].sum()-area[0:4].sum())/area[0:4].sum())

# Contribution of RF to total runoff
print(all_gt / np.sum(coeffs['runoff_gt']))











               

               

               

               

