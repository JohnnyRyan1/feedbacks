#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute empirical relationships between W m-2 and runoff in different elevation bins.

"""

# Import modules
import xarray as xr
import numpy as np
import glob
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import pandas as pd

#%%

# Define user
user = 'jryan4'

# Import data
alt_path = '/Users/jryan4/Dropbox (University of Oregon)/published/snowfall/data/'
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Define MAR files
mar_files = sorted(glob.glob(alt_path + 'mar-v3-12-1/*'))

# Derive MAR surface heights
mar = xr.open_dataset(mar_files[0])
sh = mar['SH'].values
msk = mar['MSK'].values

# Define elevations
elevations = np.arange(0, 3600, 200)
#%%

watts_by_elev = []
runoff_by_elev = []
area_by_elev = []
for e in range(len(elevations) - 1):
    print(e)
    elevation_mask = (sh > elevations[e]) & (sh < elevations[e+1]) & (msk > 90)
    
    watts_by_year = []
    runoff_by_year = []
    for m in mar_files:
        mar = xr.open_dataset(m)
        watts = mar['LWD'] + (mar['SWD']*(1-mar['AL2'][:,0,:,:])) + mar['SHF'] + mar['LHF']
        runoff = mar['RU'][:,0,:,:]      
                 
        watts_by_year.append(np.nanmean(np.nanmean(watts[152:244, :, :], axis=0)[elevation_mask]))
        runoff_by_year.append(np.nanmean(np.nanmean(runoff[152:244, :, :], axis=0)[elevation_mask]))
    
    watts_by_elev.append(watts_by_year)
    runoff_by_elev.append(runoff_by_year)
    area_by_elev.append(np.sum(elevation_mask)*10)
        
#%%

# Save coefficients as csv
coeffs = []
for i in range(len(elevations) - 1):
    coeffs.append(stats.linregress(watts_by_elev[i], runoff_by_elev[i]))

watts_df = pd.DataFrame(watts_by_elev)
runoff_df = pd.DataFrame(runoff_by_elev)

watts_by_elev = np.array(watts_by_elev)
runoff_by_elev = np.array(runoff_by_elev)
df1 = pd.DataFrame(runoff_by_elev/watts_by_elev)

df = pd.DataFrame(coeffs)
df['elev'] = elevations[:-1]
df['watts'] = watts_df.mean(axis=1)
df['runoff'] = runoff_df.mean(axis=1)

# Save as csv
df.to_csv(path + 'watts-to-runoff-coeffs.csv')
df1.to_csv(path + 'runoff-factors.csv', index=False)

#%%


# Import data
df = pd.read_csv(path + 'watts-to-runoff-coeffs.csv')

xp = np.linspace(300, 360, 100)

b = [2, 4, 8, 14]

stats_list1 = stats.linregress(watts_by_elev[b[0]], runoff_by_elev[b[0]])
stats_list2 = stats.linregress(watts_by_elev[b[1]], runoff_by_elev[b[1]])
stats_list3 = stats.linregress(watts_by_elev[b[2]], runoff_by_elev[b[2]])
stats_list4 = stats.linregress(watts_by_elev[b[3]], runoff_by_elev[b[3]])

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(11, 9))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(watts_by_elev[b[0]], runoff_by_elev[b[0]], 
            color=c2, zorder=2, s=100, alpha=0.8, label='')
ax1.plot(np.array(watts_by_elev[b[0]]), np.array(watts_by_elev[b[0]])*stats_list1[0] + stats_list1[1], '-', color='k', lw=2)
ax1.set_title('%s to %s m a.s.l.' % (elevations[b[0]], elevations[b[0]] + 200), fontsize=16)

ax2.scatter(watts_by_elev[b[1]], runoff_by_elev[b[1]], 
            color=c2, zorder=2, s=100, alpha=0.8, label='')
ax2.plot(watts_by_elev[b[1]], np.array(watts_by_elev[b[1]])*stats_list2[0] + stats_list2[1], '-', color='k', lw=2)
ax2.set_title('%s to %s m a.s.l.' % (elevations[b[1]], elevations[b[1]] + 200), fontsize=16)

ax3.scatter(watts_by_elev[b[2]], runoff_by_elev[b[2]], 
            color=c2, zorder=2, s=100, alpha=0.8, label='')
ax3.plot(watts_by_elev[b[2]], np.array(watts_by_elev[8])*stats_list3[0] + stats_list3[1], '-', color='k', lw=2)
ax3.set_title('%s to %s m a.s.l.' % (elevations[b[2]], elevations[b[2]] + 200), fontsize=16)

ax4.scatter(watts_by_elev[b[3]], runoff_by_elev[b[3]], 
            color=c2, zorder=2, s=100, alpha=0.8, label='')
#ax4.plot(watts_by_elev[12], np.array(watts_by_elev[12])*stats_list4[0] + stats_list4[1], '-', color='k', lw=2)
ax4.set_title('%s to %s m a.s.l.' % (elevations[b[3]], elevations[b[3]] + 200), fontsize=16)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list1[2]**2),
    r'Std. err. = %.2f' % stats_list1[4],
    r'Slope = %.2f mm W$^{-1}$' % stats_list1[0]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax1.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list2[2]**2),
    r'Std. err. = %.2f' % stats_list2[4],
    r'Slope = %.2f mm W$^{-1}$' % stats_list2[0]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax2.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list3[2]**2),
    r'Std. err. = %.2f' % stats_list3[4],
    r'Slope = %.2f mm W$^{-1}$' % stats_list3[0]))
text_box = AnchoredText(textstr, frameon=True, loc=1, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax3.add_artist(text_box)

# Add stats
textstr = '\n'.join((
    r'R$^{2}$ = %.2f' % (stats_list4[2]**2),
    r'Std. err. = %.2f' % stats_list4[4],
    r'Slope = %.2f mm W$^{-1}$' % stats_list4[0]))
text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5, prop=dict(size=14))
text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
plt.setp(text_box.patch, facecolor='white', alpha=0.7)
ax4.add_artist(text_box)

ax1.set_ylabel('Runoff (mm w.e.)', fontsize=14)
ax3.set_ylabel('Runoff (mm w.e.)', fontsize=14)

ax3.set_xlabel('Energy available for melt (W)', fontsize=14)
ax4.set_xlabel('Energy available for meltwater (W)', fontsize=14)

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=14)


#%%


# Temperature vs. ice sheet runoff

temp_by_year = []
runoff_by_year = []
for m in mar_files:
    mar = xr.open_dataset(m)
    runoff = mar['RU'][:,0,:,:]
    temp = mar['TT'][:,0,:,:]
    
    runoff_by_year.append(np.nansum(np.nansum(runoff[152:244, :, :], axis=0)[(msk > 90)]))
    temp_by_year.append(np.nanmean(np.nanmean(temp[152:244, :, :], axis=0)[(msk > 90)]))


runoff_by_year = np.array(runoff_by_year)
temp_by_year = np.array(temp_by_year)

# Convert from mm to Gt
runoff_gt = runoff_by_year / 1e+6 * 10 * 10
temp_anom = temp_by_year - np.mean(temp_by_year)


#%%
stats_list1 = stats.linregress(temp_anom, runoff_gt)

# Plot
fig, (ax1) = plt.subplots(nrows=1, ncols=1, figsize=(7, 7))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.scatter(temp_anom, runoff_gt, 
            color=c2, zorder=2, s=100, alpha=0.8, label='')
ax1.plot(np.sort(temp_anom), np.sort(temp_anom)*stats_list1[0] + stats_list1[1], 
         color='k', lw=2, ls='dashed')

ax1.set_xlabel('Air temperature anomaly (K)', fontsize=14)

ax1.set_ylabel('Meltwater runoff (Gt)', fontsize=14)

ax1.grid(linestyle='dotted', lw=1, zorder=1)
ax1.tick_params(axis='both', which='major', labelsize=13)
    
    
    
    
    
    
    
    
    
    