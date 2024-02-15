#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analysis at the grid cell scale

"""

# Import modules
import xarray as xr
import numpy as np
from scipy import stats
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

#%%

# Bin forcing into elevations
elevations = np.arange(0, 3400, 200)

ice, snowline, snow = [], [], []
area = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    ice.append(np.nanmean(data['ice'].values[elevation_mask]))
    snowline.append(np.nanmean(data['snowline'].values[elevation_mask]))
    snow.append(np.nanmean(data['snow'].values[elevation_mask]))

#%%

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(13, 6))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(ice, elevations[:-1], color=c1, zorder=2, alpha=0.8, label='Glacier ice')
ax1.plot(snowline, elevations[:-1], color=c2, zorder=2, alpha=0.8, label='Snowline')
ax1.plot(snow, elevations[:-1], color=c3, zorder=2, alpha=0.8, label='Snow')
ax1.set_ylim(0, 3200)
ax1.legend()

ax3.plot(ice*area, elevations[:-1], color=c1, zorder=2, alpha=0.8, label='Glacier ice')
ax3.plot(snowline*area, elevations[:-1], color=c2, zorder=2, alpha=0.8, label='Snowline')
ax3.plot(snow*area, elevations[:-1], color=c3, zorder=2, alpha=0.8, label='Snow')
ax3.set_ylim(0, 3200)
ax3.legend()
ax3.set_yticklabels([])
ax3.set_yticks([])

ax2.barh(range(len(area)), area, align='edge',  alpha=0.2, color='blue', edgecolor='k')
ax2.set_ylim(0,17)
ax2.set_yticklabels([])
ax2.set_yticks([])

#%%

"""
W vs. melt by elevation?

"""




               

               

               

               

