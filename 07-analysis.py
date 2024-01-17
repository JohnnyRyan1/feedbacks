#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analysis

"""

# Import modules
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from sklearn.metrics import mean_squared_error


#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

#%%

# Imoport data
data = pd.read_csv(path + 'results.csv')
clim = pd.read_csv(path + 'ice-sheet-climate.csv')

#%%
# Compute radiative forcing in W m-2

# SWnet due to reducing glacier ice albedo 
data['ice_feedback'] = data['observed_albedo'] - data['fixed_ice']

# SWnet due to reducing snow albedo
data['snow_feedback'] = data['observed_albedo'] - data['fixed_snow']

# SWnet due to snowline fluctuations
data['snowline_feedback'] = data['fixed_snow_ice'] - data['fixed_snow_all']

# Compute temperature anomaly
clim['t2m_anom'] = clim['t2m'] - np.mean(clim['t2m'])

#%%

# Define uncertainty scenarios
bounds = np.arange(0, data.shape[0]+22, 22)

snowline_feedback = []
snow_feedback = []
ice_feedback = []
for i in range(len(bounds) - 1):
    snowline_feedback.append(data['snowline_feedback'].iloc[bounds[i]:bounds[i+1]].values)
    snow_feedback.append(data['snow_feedback'].iloc[bounds[i]:bounds[i+1]].values)
    ice_feedback.append(data['ice_feedback'].iloc[bounds[i]:bounds[i+1]].values)


#%%

# Define years
n = np.arange(2002, 2024, 1)
xp = np.linspace(-2, 2, 100)

"""
Scatter plot showing radiative forcing in vs. temperature

"""

stats_list1 = stats.linregress(clim['t2m_anom'], snowline_feedback[4])
stats_list2 = stats.linregress(clim['t2m_anom'], snow_feedback[4])
stats_list3 = stats.linregress(clim['t2m_anom'], ice_feedback[4])

z1 = np.polyfit(clim['t2m_anom'], snowline_feedback[4], 1)
p1 = np.poly1d(z1)

z2 = np.polyfit(clim['t2m_anom'], snowline_feedback[4], 2)
p2 = np.poly1d(z2)

rmse1 = mean_squared_error(snowline_feedback[4], p1(clim['t2m_anom']), squared=False)
rmse2 = mean_squared_error(snowline_feedback[4], p2(clim['t2m_anom']), squared=False)

z3 = np.polyfit(clim['t2m_anom'], ice_feedback[4], 1)
p3 = np.poly1d(z3)

z4 = np.polyfit(clim['t2m_anom'], ice_feedback[4], 2)
p4 = np.poly1d(z4)

rmse3 = mean_squared_error(ice_feedback[4], p3(clim['t2m_anom']), squared=False)
rmse4 = mean_squared_error(ice_feedback[4], p4(clim['t2m_anom']), squared=False)

fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(clim['t2m_anom'], snowline_feedback[4], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snowline')
ax1.plot(xp, p1(xp), '-', color='k', lw=2)
ax1.plot(xp, p2(xp), '-', color='blue', lw=2)
ax1.set_title('Snowline', fontsize=16)
ax2.scatter(clim['t2m_anom'], snow_feedback[4], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snow')
ax2.set_title('Snow', fontsize=16)
ax3.scatter(clim['t2m_anom'], ice_feedback[4], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Glacier ice')
ax3.plot(xp, p3(xp), '-', color='k', lw=2)
ax3.plot(xp, p4(xp), '-', color='blue', lw=2)
ax3.set_title('Glacier ice', fontsize=16)

# =============================================================================
# for i, txt in enumerate(n):
#     ax1.annotate(txt, (clim['t2m'].values[i], snowline_feedback[4][i]), fontsize=14)
# =============================================================================

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list1[2]**2, ),
    r'RMSE = %.2f' % rmse2))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list2[2]**2, ),))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list3[2]**2, ),
    r'RMSE = %.2f' % rmse4))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)

ax1.set_ylabel('Radiative forcing (W m$^{-2}$)', fontsize=14)

for ax in [ax1, ax2, ax3]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('Summer T$_{500hPa}$ anomaly (K)', fontsize=14)
    ax.set_xlim(-1.7, 1.7)

fig.tight_layout()


#%%

"""
Scatter plot showing radiative forcing in vs. cloudiness

"""

stats_list1 = stats.linregress(clim['cloudiness'], snowline_feedback[4])
stats_list2 = stats.linregress(clim['cloudiness'], snow_feedback[4])
stats_list3 = stats.linregress(clim['cloudiness'], ice_feedback[4])


fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(clim['cloudiness'], snowline_feedback[4], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snowline')
ax1.set_title('Snowline', fontsize=16)
ax2.scatter(clim['cloudiness'], snow_feedback[4], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Snow')
ax2.set_title('Snow', fontsize=16)
ax3.scatter(clim['cloudiness'], ice_feedback[4], 
            color=c2, zorder=2, s=100, alpha=0.8, label='Glacier ice')
ax3.set_title('Glacier ice', fontsize=16)

# =============================================================================
# for i, txt in enumerate(n):
#     ax1.annotate(txt, (clim['t2m'].values[i], snowline_feedback[4][i]), fontsize=14)
# =============================================================================

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list1[2]**2, ),))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list2[2]**2, ),))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

# Add stats 
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list3[2]**2, ),))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)

ax1.set_ylabel('Radiative forcing (W m$^{-2}$)', fontsize=14)

for ax in [ax1, ax2, ax3]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlabel('Cloudiness (%)', fontsize=14)

fig.tight_layout()

#%%





