#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DESCRIPTION

Draft figures for paper

"""

# Import modules
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import colors
import matplotlib.ticker as tkr

#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

#%%
# Import data
s_forcing = xr.open_dataset(path + 'final-surface-forcing-grids.nc')
feedback = xr.open_dataset(path + 'feedback-stats.nc')
c_forcing_sw = xr.open_dataset(path + 'final-cloud-forcing-grids.nc')
c_forcing_lw = xr.open_dataset(path + 'allwave-t2m-downscaled.nc')

ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')
mask = ismip_1km['GIMP'].values

lons = s_forcing['longitude']
lats = s_forcing['latitude']

#%%

# Choose resolution
r = 10

# Resize for more convenient plotting
lons = lons[::r,::r]
lats = lats[::r,::r]

#%%
# Snowline radiative forcing
snowline_data = np.nanmean(s_forcing['snowline'], axis=2)
snowline_data[mask == 0] = np.nan

# Resize for more convenient plotting
snowline_data = snowline_data[::r,::r]

#%%
fig = plt.figure(figsize=(4, 4))
v = np.arange(0, 75, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, snowline_data, v, transform=ccrs.PlateCarree(), vmin=0, vmax=75,
             cmap='Reds', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(0, 80, 10))
cbar.ax.set_yticklabels(np.arange(0, 80, 10)) 
cbar.set_label('Radiative forcing  (W m$^{-2}$)', rotation=270, labelpad=12)
plt.savefig(savepath + 'snowline.png', dpi=300)


#%%

ice_data = np.nanmean(s_forcing['ice'], axis=2)
ice_data[mask == 0] = np.nan

# Resize for more convenient plotting
ice_data = ice_data[::r,::r]

#%%
fig = plt.figure(figsize=(4, 4))
v = np.arange(0, 75, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, ice_data, v, transform=ccrs.PlateCarree(), vmin=0, vmax=75,
             cmap='Reds', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(0, 80, 10))
cbar.ax.set_yticklabels(np.arange(0, 80, 10)) 
cbar.set_label('Radiative forcing (W m$^{-2}$)', rotation=270, labelpad=12)
plt.savefig(savepath + 'ice.png', dpi=300)

#%%
snow_data = np.nanmean(s_forcing['snow'], axis=2)
snow_data[mask == 0] = np.nan

# Resize for more convenient plotting
snow_data = snow_data[::r,::r]

#%%
fig = plt.figure(figsize=(4, 4))
v = np.arange(0, 75, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, snow_data, v, transform=ccrs.PlateCarree(), vmin=0, vmax=75,
             cmap='Reds', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(0, 80, 10))
cbar.ax.set_yticklabels(np.arange(0, 80, 10)) 
cbar.set_label('Radiative forcing  (W m$^{-2}$)', rotation=270, labelpad=12)
plt.savefig(savepath + 'snow.png', dpi=300)

#%%

# Cloud radiative forcing
sw_data = np.nanmean((c_forcing_sw['swnet_cloud'] - c_forcing_sw['swnet_no_cloud']), axis=2)
lw_data = np.nanmean((c_forcing_lw['lwd_allsky'] - c_forcing_lw['lwd_clrsky']), axis=2)
net_cloud = sw_data + lw_data
net_cloud[mask == 0] = np.nan

# Resize for more convenient plotting
net_cloud = net_cloud[::r,::r]

#%%
divnorm=colors.TwoSlopeNorm(vmin=-35, vcenter=0, vmax=75)

fig = plt.figure(figsize=(4, 4))
v = np.arange(-30, 75, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, net_cloud, v, transform=ccrs.PlateCarree(), vmin=-30, vmax=75,
             cmap='coolwarm', zorder=3, norm=divnorm)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(-30, 80, 10))
cbar.ax.set_yticklabels(np.arange(-30, 80, 10)) 
cbar.set_label('Radiative forcing  (W m$^{-2}$)', rotation=270, labelpad=12)
plt.savefig(savepath + 'clouds.png', dpi=300)


#%%

# Feedback strength
surface = feedback['surface_feedback'].values[:,:,0]

surface[mask == 0] = np.nan

# Resize for more convenient plotting
surface = surface[::r,::r]

divnorm=colors.TwoSlopeNorm(vmin=-5., vcenter=0., vmax=25)

fig = plt.figure(figsize=(4, 4))
v = np.arange(-5, 25, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, surface, v, transform=ccrs.PlateCarree(), vmin=-5, vmax=25,
             cmap='coolwarm', zorder=3, norm=divnorm)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(-5, 25, 5))
cbar.ax.set_yticklabels(np.arange(-5, 25, 5)) 
cbar.set_label('Radiative feedback  (W m$^{-2}$ K$^{-1}$)', rotation=270, labelpad=12)
plt.savefig(savepath + 'surface-feedback-slope.png', dpi=300)

#%%

# Feedback significances
p = feedback['surface_feedback'].values[:,:,3] > 0.1
p = p.astype('float')
p[mask == 0] = np.nan

# Resize for more convenient plotting
p = p[::r,::r]

fig = plt.figure(figsize=(4, 4))
#v = np.arange(-5, 25, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, p, transform=ccrs.PlateCarree(),
             cmap='Blues_r', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
plt.savefig(savepath + 'surface-feedback-sig.png', dpi=300)


#%%

# Feedback significances
p = feedback['cloud_net_feedback'].values[:,:,3] > 0.1
p = p.astype('float')
p[mask == 0] = np.nan

# Resize for more convenient plotting
p = p[::r,::r]

fig = plt.figure(figsize=(4, 4))
#v = np.arange(-5, 25, 1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, p, transform=ccrs.PlateCarree(),
             cmap='Blues_r', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
plt.savefig(savepath + 'cloud-feedback-sig.png', dpi=300)


#%%

t2m_r = np.nanmean(c_forcing_lw['t2m_r'], axis=2)

t2m_r[mask == 0] = np.nan

# Resize for more convenient plotting
t2m_r = t2m_r[::r,::r]

fig = plt.figure(figsize=(4, 4))
v = np.arange(-1, 0.1, 0.05)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, t2m_r, v, transform=ccrs.PlateCarree(), vmin=-1, vmax=0,
             cmap='Blues_r', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(-1, 0.1, 0.2))
cbar.ax.set_yticklabels(np.arange(-1, 0.1, 0.2)) 
cbar.ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%.1f'))
cbar.set_label('Correlation with elevation', rotation=270, labelpad=12)
plt.savefig(savepath + 'downscale-t2m-correlations.png', dpi=300)

#%%
swd_allsky_r = np.nanmean(c_forcing_lw['swd_allsky_r'], axis=2)

swd_allsky_r[mask == 0] = np.nan

# Resize for more convenient plotting
swd_allsky_r = swd_allsky_r[::r,::r]

divnorm=colors.TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)

fig = plt.figure(figsize=(4, 4))
v = np.arange(-1, 1.1, 0.1)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, swd_allsky_r, v, transform=ccrs.PlateCarree(), vmin=-1, vmax=1,
             cmap='coolwarm', zorder=3, norm=divnorm)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(-1, 1.2, 0.4))
cbar.ax.set_yticklabels(np.arange(-1, 1.2, 0.4)) 
cbar.ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%.1f'))
cbar.set_label('Correlation with elevation', rotation=270, labelpad=12)
plt.savefig(savepath + 'downscale-swd-allsky-correlations.png', dpi=300)


#%%
lwd_allsky_r = np.nanmean(c_forcing_lw['lwd_allsky_r'], axis=2)

lwd_allsky_r[mask == 0] = np.nan

# Resize for more convenient plotting
lwd_allsky_r = lwd_allsky_r[::r,::r]

fig = plt.figure(figsize=(4, 4))
v = np.arange(-1, 0.1, 0.05)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
plt.contourf(lons, lats, lwd_allsky_r, v, transform=ccrs.PlateCarree(), vmin=-1, vmax=0,
             cmap='Blues_r', zorder=3)
ax.coastlines(resolution='50m', color='black', linewidth=0.5)
ax.add_feature(cfeature.LAND, facecolor='antiquewhite', zorder=2)
ax.add_feature(cfeature.OCEAN,facecolor='white', zorder=1)
cbar = plt.colorbar(ticks=np.arange(-1, 0.1, 0.2))
cbar.ax.set_yticklabels(np.arange(-1, 0.1, 0.2)) 
cbar.ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%.1f'))
cbar.set_label('Correlation with elevation', rotation=270, labelpad=12)
plt.savefig(savepath + 'downscale-lwd-allsky-correlations.png', dpi=300)

#%%















