#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plot MAR values

"""

# Import modules
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import pandas as pd

#%%
# Define user
user = 'jryan4'

# Import data
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

# Read W to mm conversion coefficients
factor = pd.read_csv(path + 'melt-factors.csv')
watts = pd.read_csv(path + 'watts.csv')
melt = pd.read_csv(path + 'melt.csv')
mar = pd.read_csv(path + 'mar-vs-temp.csv')
coeffs = pd.read_csv(path + 'final-coeffs.csv')

# Define elevations
elevations = np.arange(0, 3600, 200)
area = np.array((9806,  12752,  18442,  27446,  42198,  58465,  74097,  93539,
       113645, 138732, 162447, 181550, 204249, 195323, 172527, 102933,
         4497))

mean_melt = melt.mean(axis=1)/1e+06*area*92
std_melt = melt.std(axis=1)/1e+06*area*92

coeffs['bulk'] = coeffs['ice_gt']+coeffs['snowline_gt']+coeffs['snow_gt']
coeffs['bulk_std'] = coeffs['ice_gt_std']+coeffs['snowline_gt_std']+coeffs['snow_gt_std']
#%%


fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

#ax1.scatter(mean_melt, elevations[:-1], color=c1, zorder=2, s=100, alpha=0.8, label='')
ax1.plot(mean_melt, elevations[:-1], color=c1, zorder=2, lw=2, alpha=0.8, 
         ls='dashed', label='Melt total')
ax1.fill_betweenx(elevations[:-1], mean_melt - std_melt, 
                  mean_melt + std_melt, zorder=1,
                  color=c1, alpha=0.2)

#ax1.scatter(coeffs['bulk'], elevations[:-1], color=c2, zorder=2, s=100, alpha=0.8, label='')
ax1.plot(coeffs['bulk'], elevations[:-1], color=c2, zorder=2, lw=2, alpha=0.8, 
         ls='dashed', label='Melt due surface radiative forcing')
ax1.fill_betweenx(elevations[:-1], coeffs['bulk'] - coeffs['bulk_std'] , 
                  coeffs['bulk'] + coeffs['bulk_std'], zorder=1,
                  color=c2, alpha=0.2)
ax1.set_ylim(0, 3200)

ax1.set_xlabel('Summer meltwater production (Gt yr$^{-1}$)', fontsize=14)
ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)
ax1.axhline(y=1600, ls='dashed', color='k', lw=2, zorder=0, alpha=0.8)
ax1.legend(fontsize=13)
ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)

fig.savefig(savepath + 'melt-vs-elevation.png', dpi=200)


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

ax1.set_xlabel('Melt factor (mm W$^{-1}$)', fontsize=14)
ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)

fig.savefig(savepath + 'watts-to-melt.png', dpi=200)

#%%

# Bin correlations into elevations
elevations = np.arange(0, 3600, 200)

correlation = []
for e in range(watts.shape[0]):
    stats_list1 = stats.linregress(watts.iloc[e], melt.iloc[e])
    correlation.append(stats_list1[2])

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 6), sharey=True)

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(correlation, elevations[:-1], color=c1, zorder=2, s=100, alpha=0.8, label='Correlation (r)')
ax1.plot(correlation, elevations[:-1], color=c1, zorder=2, lw=2, alpha=0.5, 
         ls='dashed', label='')
ax1.set_ylim(0, 3400)

ax2.scatter(melt.mean(axis=1)*92/1000, elevations[:-1], color=c2, zorder=2, s=100, 
            alpha=0.8, label='Melt (mm)')
ax2.plot(melt.mean(axis=1)*92/1000, elevations[:-1], color=c2, zorder=2, lw=2, alpha=0.5, 
         ls='dashed', label='')

ax1.set_xlabel('Correlation coefficient (r)', fontsize=14)
ax2.set_xlabel('Summer melt (m)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)
ax2.grid(linestyle='dotted', lw=1, zorder=1)
ax2.tick_params(axis='both', which='major', labelsize=13)
ax1.text(0.90, 0.93, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.90, 0.91, "b", fontsize=24, transform=ax2.transAxes)


fig.savefig(savepath + 'coefficients-vs-elevation.png', dpi=200)


#%%

stats_list1 = stats.linregress(mar['temp'], mar['melt'])
xp = np.linspace(np.sort(mar['temp']).min() - 1, np.sort(mar['temp']).max() + 1, mar['temp'].shape[0])

# Plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(mar['temp'], mar['melt'], 
            color=c2, zorder=2, s=100, alpha=0.8, label='')
ax1.plot(xp, xp*stats_list1[0] + stats_list1[1], 
         color='k', lw=2, ls='dashed')

ax1.set_xlabel('Air temperature (2 m) anomaly (K)', fontsize=14)
ax1.set_ylabel('Meltwater production (Gt)', fontsize=14)
ax1.set_xlim(-1.2, 1.7)
ax1.set_ylim(260, 775)
ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)
ax1.axhline(y=396, c='k', ls='dotted', lw=2)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list1[2]**2),
    r'Std. err. = %.2f Gt' % stats_list1[4]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)


fig.savefig(savepath + 'temp-vs-melt.png', dpi=200)

#%%






