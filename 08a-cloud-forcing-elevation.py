#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Cloud radiative forcing and feedback analysis.

"""

# Import modules
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'
alt_path = '/Users/' + user + '/Dropbox (University of Oregon)/published/clouds/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

# Read W to mm conversion coefficients
coeffs = pd.read_csv(path + 'watts-to-melt-coeffs.csv')

# Import data
data = xr.open_dataset(path + 'feedback-stats.nc')
cloud_sw = xr.open_dataset(path + 'final-cloud-forcing-grids.nc')
swd_lwd = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')

#%%

xr_stats_sw_sig = data['cloud_sw_feedack'].copy()
xr_stats_sw_sig = xr_stats_sw_sig.where(xr_stats_sw_sig[:,:,3] < 0.1)
c_mask_sw = np.isfinite(xr_stats_sw_sig[:,:,0]).values
c_mask_sw[mask == False] = False

xr_stats_lw_sig = data['cloud_lw_feedback'].copy()
xr_stats_lw_sig = xr_stats_lw_sig.where(xr_stats_lw_sig[:,:,3] < 0.1)
c_mask_lw = np.isfinite(xr_stats_lw_sig[:,:,0]).values
c_mask_lw[mask == False] = False

xr_stats_net_sig = data['cloud_net_feedback'].copy()
xr_stats_net_sig = xr_stats_net_sig.where(xr_stats_net_sig[:,:,3] < 0.1)
c_mask_net = np.isfinite(xr_stats_net_sig[:,:,0]).values
c_mask_net[mask == False] = False

#%%

# Bin forcing into elevations
elevations = np.arange(0, 3600, 200)

clouds_sw, clouds_lw = [], []
clouds_sw_std, clouds_lw_std = [], []
area = []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area.append(elevation_mask.sum())
    
    cre_sw_mean = np.nanmean((cloud_sw['swnet_no_cloud'].values - cloud_sw['swnet_cloud'].values), axis=2)
    cre_sw_std = np.nanstd((cloud_sw['swnet_no_cloud'].values - cloud_sw['swnet_cloud'].values), axis=2)
    
    cre_lw_mean = np.nanmean((swd_lwd['lwd_allsky'].values - swd_lwd['lwd_clrsky'].values), axis=2)
    cre_lw_std = np.nanstd((swd_lwd['lwd_allsky'].values - swd_lwd['lwd_clrsky'].values), axis=2)
    
    clouds_sw.append(np.nanmean(cre_sw_mean[elevation_mask]))
    clouds_sw_std.append(np.nanmean(cre_sw_std[elevation_mask]))
    clouds_lw.append(np.nanmean(cre_lw_mean[elevation_mask]))
    clouds_lw_std.append(np.nanmean(cre_lw_std[elevation_mask]))
                         
area = np.array(area)

clouds_sw = np.array(clouds_sw)
clouds_lw = np.array(clouds_lw)
clouds_sw_std = np.array(clouds_sw_std)
clouds_lw_std = np.array(clouds_lw_std)

#%%

# Compute melt factor
coeffs['factor'] = coeffs['melt'] / coeffs['watts']

coeffs['clouds_sw'] = coeffs['factor'] * clouds_sw
coeffs['clouds_lw'] = coeffs['factor'] * clouds_lw

coeffs['clouds_sw_std'] = coeffs['factor'] * clouds_sw_std
coeffs['clouds_lw_std'] = coeffs['factor'] * clouds_lw_std

coeffs['cloud_sw_gt'] = coeffs['clouds_sw'] / 1e+06 * area * 92
coeffs['cloud_lw_gt'] = coeffs['clouds_lw'] / 1e+06 * area * 92

coeffs['cloud_sw_gt_std'] = coeffs['clouds_sw_std'] / 1e+06 * area * 92
coeffs['cloud_lw_gt_std'] = coeffs['clouds_lw_std'] / 1e+06 * area * 92

#%%

# Bin feedbacks into elevations
elevations = np.arange(0, 3600, 200)

area_all, area_valid_sw, area_valid_lw, area_valid_net = [], [], [], []
bulk_sw, bulk_lw, bulk_net = [], [], []
bulk_sw_std, bulk_lw_std, bulk_net_std = [], [], []

for e in range(len(elevations) - 1):
    elevation_mask = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1])
    area_all.append(elevation_mask.sum())
  
    cloud_feedback = data['cloud_sw_feedack'].values[:,:,0]
    cloud_feedback[c_mask_sw == False] = 0
    bulk_sw.append(np.nanmean(cloud_feedback[elevation_mask]))
    bulk_sw_std.append(np.nanstd(cloud_feedback[elevation_mask]))
    
    cloud_feedback = data['cloud_lw_feedback'].values[:,:,0]
    cloud_feedback[c_mask_lw == False] = 0
    bulk_lw.append(np.nanmean(cloud_feedback[elevation_mask]))
    bulk_lw_std.append(np.nanstd(cloud_feedback[elevation_mask]))
    
    cloud_feedback = data['cloud_net_feedback'].values[:,:,0]
    cloud_feedback[c_mask_lw == False] = 0
    bulk_net.append(np.nanmean(cloud_feedback[elevation_mask]))
    bulk_net_std.append(np.nanstd(cloud_feedback[elevation_mask]))

    area_valid_c = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (c_mask_sw == True)
    area_valid_sw.append(area_valid_c.sum())
    
    area_valid_c = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (c_mask_lw == True)
    area_valid_lw.append(area_valid_c.sum())   
    
    area_valid_n = (mask == True) & (ismip_1km['SRF'].values > elevations[e]) & (ismip_1km['SRF'].values < elevations[e+1]) & (c_mask_net == True)
    area_valid_net.append(area_valid_n.sum())   

area_all = np.array(area_all)
bulk_sw = np.array(bulk_sw)
bulk_lw = np.array(bulk_lw)
bulk_net = np.array(bulk_net)

coeffs['bulk_sw'] = coeffs['factor'] * bulk_sw
coeffs['bulk_sw_std'] = coeffs['factor'] * bulk_sw_std
coeffs['bulk_sw_gt'] = coeffs['bulk_sw'] / 1e+06 * area * 92
coeffs['bulk_sw_gt_std'] = coeffs['bulk_sw_std'] / 1e+06 * area * 92

coeffs['bulk_lw'] = coeffs['factor'] * bulk_lw
coeffs['bulk_lw_std'] = coeffs['factor'] * bulk_lw_std
coeffs['bulk_lw_gt'] = coeffs['bulk_lw'] / 1e+06 * area * 92
coeffs['bulk_lw_gt_std'] = coeffs['bulk_lw_std'] / 1e+06 * area * 92

coeffs['bulk_net'] = coeffs['factor'] * bulk_net
coeffs['bulk_net_std'] = coeffs['factor'] * bulk_lw_std
coeffs['bulk_net_gt'] = coeffs['bulk_net'] / 1e+06 * area * 92
coeffs['bulk_net_gt_std'] = coeffs['bulk_net_std'] / 1e+06 * area * 92
#%%

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

# Define colour map
c1 = '#E05861'
c2 = '#616E96'
c3 = '#F8A557'
c4 = '#3CBEDD'

ax1.plot(-clouds_sw, elevations[:-1], color=c2, zorder=2, lw=2, 
         alpha=0.8, label='SW')
ax1.fill_betweenx(elevations[:-1],
                 -clouds_sw+clouds_sw_std, -clouds_sw-clouds_sw_std,
                 zorder=1, color=c2, alpha=0.2)
ax1.plot(clouds_lw, elevations[:-1], color=c1, zorder=2, lw=2, 
         alpha=0.8, label='LW')
ax1.fill_betweenx(elevations[:-1],
                 clouds_lw+clouds_lw_std, clouds_lw-clouds_lw_std,
                 zorder=1, color=c1, alpha=0.2)
ax1.plot(clouds_lw-clouds_sw, elevations[:-1], color=c4, zorder=2, lw=2, 
         alpha=0.8, label='Net')
ax1.fill_betweenx(elevations[:-1],
                 clouds_lw-clouds_sw+np.sqrt(clouds_lw_std**2 + clouds_sw_std**2), 
                 clouds_lw-clouds_sw-np.sqrt(clouds_lw_std**2 + clouds_sw_std**2),
                 zorder=1, color=c4, alpha=0.2)

ax1.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax1.axvline(x=0, ls='dashed', color='k', zorder=1, alpha=0.5)

ax1.legend(fontsize=13, loc=1)
ax1.set_xlim(-40, 40)
ax1.set_ylim(0, 3400)

ax2.plot(coeffs['cloud_lw_gt'] - coeffs['cloud_sw_gt'], elevations[:-1], 
         color=c4, zorder=2, alpha=0.8, label='')
ax2.fill_betweenx(elevations[:-1],
                 coeffs['cloud_lw_gt'] - coeffs['cloud_sw_gt']+
                 np.sqrt(coeffs['cloud_lw_gt_std']**2 + coeffs['cloud_sw_gt_std']**2), 
                 coeffs['cloud_lw_gt'] - coeffs['cloud_sw_gt']-
                 np.sqrt(coeffs['cloud_lw_gt_std']**2 + coeffs['cloud_sw_gt_std']**2),
                 zorder=1, color=c4, alpha=0.2)

ax2.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax2.set_ylim(0, 3400)
ax2.axvline(x=0, ls='dashed', color='k', zorder=1, alpha=0.5)
#ax2.legend(fontsize=13)
ax2.set_yticklabels([])

ax3.plot(bulk_sw, elevations[:-1], color=c2, zorder=2, alpha=0.8, label='SW')
ax3.fill_betweenx(elevations[:-1],
                 bulk_sw+bulk_sw_std, bulk_sw-bulk_sw_std,
                 zorder=1, color=c2, alpha=0.2)
ax3.plot(bulk_lw, elevations[:-1], color=c1, zorder=2, alpha=0.8, label='LW')
ax3.fill_betweenx(elevations[:-1],
                 bulk_lw+bulk_lw_std, bulk_lw-bulk_lw_std,
                 zorder=1, color=c1, alpha=0.2)
ax3.plot(bulk_sw+bulk_lw, elevations[:-1], color=c4, zorder=2, alpha=0.8, label='Net')
ax3.fill_betweenx(elevations[:-1],
                 bulk_net+bulk_net_std, bulk_net-bulk_net_std,
                 zorder=1, color=c4, alpha=0.2)

ax3.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.axvline(x=0, ls='dashed', color='k', zorder=1, alpha=0.5)
ax3.set_ylim(0, 3400)
ax3.legend(fontsize=13)
ax3.legend(fontsize=13, loc=1)


ax4.plot(coeffs['bulk_sw_gt'], elevations[:-1], color=c2, zorder=2, 
         lw=2, alpha=0.8, label='SW')
ax4.fill_betweenx(elevations[:-1],
                 coeffs['bulk_sw_gt']+coeffs['bulk_sw_gt_std'], coeffs['bulk_sw_gt']-coeffs['bulk_sw_gt_std'],
                 zorder=1, color=c2, alpha=0.2)
ax4.plot(coeffs['bulk_lw_gt'], elevations[:-1], color=c1, 
         lw=2, zorder=2, alpha=0.8, label='LW')
ax4.fill_betweenx(elevations[:-1],
                 coeffs['bulk_lw_gt']+coeffs['bulk_lw_gt_std'], coeffs['bulk_lw_gt']-coeffs['bulk_lw_gt_std'],
                 zorder=1, color=c1, alpha=0.2)
ax4.plot(coeffs['bulk_sw_gt'] + coeffs['bulk_lw_gt'], elevations[:-1], color=c4, 
         lw=2, zorder=2, alpha=0.8, label='Net')
ax4.fill_betweenx(elevations[:-1],
                 coeffs['bulk_net_gt']+coeffs['bulk_net_gt_std'], coeffs['bulk_net_gt']-coeffs['bulk_net_gt_std'],
                 zorder=1, color=c4, alpha=0.2)
ax4.axvline(x=0, ls='dashed', color='k', zorder=1, alpha=0.5)
ax4.axhline(y=1600, ls='dashed', color='k', zorder=1, alpha=0.5)
ax4.set_ylim(0, 3400)
ax4.legend(fontsize=13)
ax4.set_yticklabels([])

ax1.set_xlabel('Radiative forcing (W m$^{-2}$)', fontsize=14)
ax2.set_xlabel('Meltwater production (Gt yr$^{-1}$)', fontsize=14)
ax3.set_xlabel('Radiative feedback (W m$^{-2}$ K$^{-1}$)', fontsize=14)
ax4.set_xlabel('Meltwater production (Gt yr$^{-1}$ K$^{-1}$)', fontsize=14)

ax1.set_ylabel('Elevation (m a.s.l.)', fontsize=14)
ax3.set_ylabel('Elevation (m a.s.l.)', fontsize=14)

for ax in [ax1, ax2, ax3, ax4]:
    ax.grid(linestyle='dotted', lw=1, zorder=1)
    ax.tick_params(axis='both', which='major', labelsize=13)

ax1.text(0.03, 0.89, "a", fontsize=24, transform=ax1.transAxes)
ax2.text(0.03, 0.87, "b", fontsize=24, transform=ax2.transAxes)
ax3.text(0.03, 0.89, "c", fontsize=24, transform=ax3.transAxes)
ax4.text(0.03, 0.89, "d", fontsize=24, transform=ax4.transAxes)

fig.savefig(savepath + 'cloud-radiative-forcing.png', dpi=200)


#%%

"""
Clouds
"""

# Decrease in SWnet due to clouds
cre_sw = np.nanmean(cloud_sw['swnet_no_cloud'] - cloud_sw['swnet_cloud'], axis=2)
print(np.nanmean(cre_sw[(mask == True)]))

cre_sw_std = np.nanstd(cloud_sw['swnet_no_cloud'] - cloud_sw['swnet_cloud'], axis=2)
print(np.nanmean(cre_sw_std[(mask == True)]))

# Increase LW due to clouds
cre_lw = np.nanmean((swd_lwd['lwd_allsky'].values - swd_lwd['lwd_clrsky'].values), axis=2)
print(np.nanmean(cre_lw[(mask == True)]))

cre_lw_std = np.nanstd((swd_lwd['lwd_allsky'].values - swd_lwd['lwd_clrsky'].values), axis=2)
print(np.nanmean(cre_lw_std[(mask == True)]))

# Cloud shading effect between 0-600 m
print(np.mean(clouds_sw[0:3]), '+/-', np.mean(clouds_sw_std[0:3]))

# Cloud shading effect >1600 m
print(np.mean(clouds_sw[8:]), '+/-', np.mean(clouds_sw_std[8:]))

# Increase in SWnet due to no clouds <1600 m a.s.l.
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 1600)
cre_sw_1600 = np.nanmean(cre_sw[elevation_mask])
cre_lw_1600 = np.nanmean(cre_lw[elevation_mask])
cre_sw_1600_std = np.nanmean(cre_sw_std[elevation_mask])
cre_lw_1600_std = np.nanmean(cre_lw_std[elevation_mask])
print(cre_lw_1600 - cre_sw_1600)
print(np.sqrt(cre_sw_1600_std**2 + cre_lw_1600_std**2))

# Reduction of meltwater melt due to clouds
print(coeffs['cloud_lw_gt'].sum() - coeffs['cloud_sw_gt'].sum())
print(np.sqrt(coeffs['cloud_lw_gt_std'].sum()**2 + coeffs['cloud_sw_gt_std'].sum()**2))
 
 #%%

# Fraction significant correlations between air temp and CRE NET for entire ice sheet
print(c_mask_net.sum()/mask.sum())

# Cloud allwave radiative feedback at the ice sheet scale
net_cloud_feedback = data['cloud_net_feedback'].values[:,:,0]
net_cloud_feedback[c_mask_net == False] = 0
print(np.nanmean(net_cloud_feedback))


#%%

# Make map
cloud_mask = c_mask_net
cloud_mask = cloud_mask.astype(int)
cloud_mask[(c_mask_sw == 1) & (cloud_mask == 1)] = 2
cloud_mask[(c_mask_sw == 1) & (cloud_mask == 0)] = 3
cloud_mask[(c_mask_lw == 1) & (cloud_mask == 1)] = 4
cloud_mask[(c_mask_lw == 1) & (cloud_mask == 0)] = 5

plt.imshow(cloud_mask, cmap='Set1')
plt.colorbar()

#%%
# Fraction of significant correlations between air temp and CRE SW for <600m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 600)
print(c_mask_sw[elevation_mask].sum()/mask[elevation_mask].sum())

# SW radiative forcing due to clouds <1600m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 600)
sw_cloud_feedback = data['cloud_sw_feedack'].values[:,:,0]
sw_cloud_feedback[c_mask_sw == False] = 0
print(np.nanmean(sw_cloud_feedback[elevation_mask]))
print(np.nanstd(sw_cloud_feedback[elevation_mask]))

# Fraction of significant correlations between air temp and CRE SW for <600m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 600)
print(c_mask_lw[elevation_mask].sum()/mask[elevation_mask].sum())

# SW radiative forcing due to clouds <1600m
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 0) & (ismip_1km['SRF'].values <= 600)
lw_cloud_feedback = data['cloud_lw_feedback'].values[:,:,0]
lw_cloud_feedback[c_mask_lw == False] = 0
print(np.nanmean(lw_cloud_feedback[elevation_mask]))
print(np.nanstd(lw_cloud_feedback[elevation_mask]))

print(np.nanmean(sw_cloud_feedback[elevation_mask]) + np.nanmean(lw_cloud_feedback[elevation_mask]))
print(np.sqrt(np.nanstd(lw_cloud_feedback[elevation_mask]**2 + sw_cloud_feedback[elevation_mask]**2)))

#%%
elevation_mask = (mask == True) & (ismip_1km['SRF'].values > 600) & (ismip_1km['SRF'].values <= 1600)
print(np.nanmean(sw_cloud_feedback[elevation_mask]))
print(np.nanmean(lw_cloud_feedback[elevation_mask]))

print(np.nanstd(sw_cloud_feedback[elevation_mask]))
print(np.nanstd(lw_cloud_feedback[elevation_mask]))


print(np.sqrt(np.nanstd(lw_cloud_feedback[elevation_mask]**2 + sw_cloud_feedback[elevation_mask]**2)))

# Meltwater due to cloud radiative feedback
print((coeffs['bulk_lw_gt'] + coeffs['bulk_sw_gt']).sum())

#%%


    






