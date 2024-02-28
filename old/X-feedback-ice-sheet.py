#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analysis at the ice sheet scale

"""

# Import modules
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

#%%

# Imoport data
pos = pd.read_csv(path + 'positive-forcing-results.csv')
neg = pd.read_csv(path + 'negatiave_forcing_results.csv')
clim = pd.read_csv(path + 'ice-sheet-climate.csv')

#%%

# Compute temperature anomaly
clim['t2m_anom'] = clim['t2m'] - np.mean(clim['t2m'])
clim['c_anom'] = (clim['cloudiness']*100) - (np.mean(clim['cloudiness'])*100)

#%%

# Define uncertainty scenarios
bounds = np.arange(0, pos.shape[0]+22, 22)

snowline_forcing = []
snow_forcing = []
ice_forcing = []
for i in range(len(bounds) - 1):
    snowline_forcing.append(pos['snowline_forcing'].iloc[bounds[i]:bounds[i+1]].values)
    snow_forcing.append(pos['snow_forcing'].iloc[bounds[i]:bounds[i+1]].values)
    ice_forcing.append(pos['ice_forcing'].iloc[bounds[i]:bounds[i+1]].values)
    

#%%

# Define years
n = np.arange(2002, 2024, 1)
xp = np.linspace(-2, 2, 100)

"""
Scatter plot showing radiative forcing vs. temperature

"""

stats_list1 = stats.linregress(clim['t2m_anom'], ice_forcing[0])
stats_list2 = stats.linregress(clim['t2m_anom'], snowline_forcing[0])
stats_list3 = stats.linregress(clim['t2m_anom'], snow_forcing[0])
stats_list4 = stats.linregress(clim['t2m_anom'], (ice_forcing[0] + snowline_forcing[0] + snow_forcing[0]))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(clim['t2m_anom'], ice_forcing[0], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Glacier ice')
ax1.plot(xp, xp*stats_list1[0] + stats_list1[1], '-', color='k', lw=2)
ax1.set_title('Glacier ice', fontsize=16)
ax2.scatter(clim['t2m_anom'], snowline_forcing[0], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snowline')
ax2.plot(xp, xp*stats_list2[0] + stats_list2[1], '-', color='k', lw=2)
ax2.set_title('Snowline', fontsize=16)
ax3.scatter(clim['t2m_anom'], snow_forcing[0], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snow')
ax3.plot(xp, xp*stats_list3[0] + stats_list3[1], '-', color='k', lw=2)
ax3.set_title('Snow', fontsize=16)
ax4.scatter(clim['t2m_anom'], (ice_forcing[0] + snowline_forcing[0] + snow_forcing[0]), 
            color=c1, zorder=2, s=100, alpha=0.8, label='All')
ax4.plot(xp, xp*stats_list4[0] + stats_list4[1], '-', color='k', lw=2)
ax4.set_title('Glacier ice + snowline + snow', fontsize=16)

# =============================================================================
# for i, txt in enumerate(n):
#     ax1.annotate(txt, (clim['t2m'].values[i], snowline_feedback[0][i]), fontsize=14)
# =============================================================================

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list1[2]**2),
    r'Std. err. = %.2f' % stats_list1[4],
    r'Slope = %.2f W m$^{-2}$ K$^{-1}$' % stats_list1[0]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list2[2]**2),
    r'Std. err. = %.2f' % stats_list2[4],
    r'Slope = %.2f W m$^{-2}$ K$^{-1}$' % stats_list2[0]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list3[2]**2),
    r'Std. err. = %.2f' % stats_list3[4],
    r'Slope = %.2f W m$^{-2}$ K$^{-1}$' % stats_list3[0]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list4[2]**2),
    r'Std. err. = %.2f' % stats_list4[4],
    r'Slope = %.2f W m$^{-2}$ K$^{-1}$' % stats_list4[0]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax4.add_artist(text_box)

ax1.set_ylabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax3.set_ylabel('Radiative forcing (W m$^{-2}$)', fontsize=14)

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim(-1.7, 1.7)

for ax in [ax3, ax4]:
    ax.set_xlabel('Summer T$_{2m}$ anomaly (K)', fontsize=14)

ax1.set_ylim(0.2, 5.2)
ax2.set_ylim(0.2, 5.2)
ax3.set_ylim(5.2, 20.2)
ax4.set_ylim(5.2, 20.2)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.89, "d", fontsize=24, transform=ax4.transAxes)

fig.tight_layout()
fig.savefig(savepath + 'radiative-forcing-temperature.png', dpi=200)

#%%

"""
Scatter plot showing radiative forcing vs. cloudiness

"""

stats_list1 = stats.linregress(clim['cloudiness'], ice_forcing[0])
stats_list2 = stats.linregress(clim['cloudiness'], snowline_forcing[0])
stats_list3 = stats.linregress(clim['cloudiness'], snow_forcing[0])
stats_list4 = stats.linregress(clim['cloudiness'], (ice_forcing[0] + snowline_forcing[0] + snow_forcing[0]))


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(clim['cloudiness'], ice_forcing[0], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snowline')
ax1.set_title('Glacier ice', fontsize=16)
ax2.scatter(clim['cloudiness'], snowline_forcing[0], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snow')
ax2.set_title('Snowline', fontsize=16)
ax3.scatter(clim['cloudiness'], snow_forcing[0], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Glacier ice')
ax3.set_title('Snow', fontsize=16)
ax4.scatter(clim['cloudiness'], (ice_forcing[0] + snowline_forcing[0] + snow_forcing[0]), 
            color=c1, zorder=2, s=100, alpha=0.8, label='All')
ax4.set_title('Glacier ice + snowline + snow', fontsize=16)

# =============================================================================
# for i, txt in enumerate(n):
#     ax1.annotate(txt, (clim['t2m'].values[i], snowline_feedback[0][i]), fontsize=14)
# =============================================================================

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list1[2]**2),
    r'Std. err. = %.2f' % stats_list1[4],
    r'Slope = %.2f W m$^{-2}$ ‰$^{-1}$' % (stats_list1[0]/10)))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list2[2]**2),
    r'Std. err. = %.2f' % stats_list2[4],
    r'Slope = %.2f W m$^{-2}$ ‰$^{-1}$' % (stats_list2[0]/10)))
text_box = AnchoredText(textstr, frameon=True, loc=3, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list3[2]**2),
    r'Std. err. = %.2f' % stats_list3[4],
    r'Slope = %.2f W m$^{-2}$ ‰$^{-1}$' % (stats_list3[0]/10)))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list4[2]**2),
    r'Std. err. = %.2f' % stats_list4[4],
    r'Slope = %.2f W m$^{-2}$ ‰$^{-1}$' % (stats_list4[0]/10)))
text_box = AnchoredText(textstr, frameon=True, loc=3, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax4.add_artist(text_box)


ax1.set_ylabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax3.set_ylabel('Radiative forcing (W m$^{-2}$)', fontsize=14)

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=14)

ax3.set_xlabel('Cloudiness anomaly (%)', fontsize=14)
ax4.set_xlabel('Cloudiness anomaly (%)', fontsize=14)

ax1.set_ylim(0.2, 6.3)
ax2.set_ylim(0.2, 6.3)
ax3.set_ylim(4.2, 14)
ax4.set_ylim(3.2, 20.2)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.89, "d", fontsize=24, transform=ax4.transAxes)


fig.tight_layout()
fig.savefig(savepath + 'radiative-forcing-clouds.png', dpi=200)


#%%

"""
Standardized linear regresion to compute relative strength of temperature vs. clouds

"""

# Standardize
t = (clim['t2m'] - np.mean(clim['t2m'])) / np.std(clim['t2m'])
c = (clim['cloudiness'] - np.mean(clim['cloudiness'])) / np.std(clim['cloudiness'])

stats1 = stats.linregress(t, (ice_forcing[0] + snowline_forcing[0] + snow_forcing[0]))
stats2 = stats.linregress(c, (ice_forcing[0] + snowline_forcing[0] + snow_forcing[0]))


































