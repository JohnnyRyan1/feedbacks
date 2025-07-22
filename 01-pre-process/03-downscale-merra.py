#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Downscale MERRA-2 longwave and shortwave downward radiation.

"""

# Import modules
import xarray as xr
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import glob
from pyresample import kd_tree, geometry
import netCDF4

#%%


# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Define mask
mask = ismip_1km['GIMP'].values

# Define target projection
target_def = geometry.SwathDefinition(lons=ismip_1km['lon'].values, lats=ismip_1km['lat'].values)

# Define MERRA-2 elevations
elev_file = xr.open_dataset(path + 'MERRA2_101.const_2d_asm_Nx.00000000.nc4.nc4')

# Define MERRA files
merra_sw_files = sorted(glob.glob(path + 'merra-swd/*.nc'))
merra_lw_files = sorted(glob.glob(path + 'merra-lwd/*.nc'))
merra_tm_files = sorted(glob.glob(path + 'merra-t2m/*.nc'))

# Convert geopotential height to elevation
elev = elev_file['PHIS'].values[0,:,:] / 9.81

# Read MERRA data to find lat/lons
merra = xr.open_dataset(merra_sw_files[0])

# Extract lat/lon
merra_lon, merra_lat = np.meshgrid(merra['longitude'].values, merra['latitude'].values)

# Define source projection
source_def = geometry.SwathDefinition(lons=merra_lon, lats=merra_lat)

#%%

def produce_stats(climatology):
        
    p, a, b, r = np.zeros((55,105)), np.zeros((55,105)), np.zeros((55,105)), np.zeros((55,105))
    
    for i in range(b.shape[0]):
        for j in range(b.shape[1]):
            
            imin = i-2
            imax = i+3
            jmin = j-2
            jmax = j+3
            
            if imin < 0:
                imin=0
            if jmin < 0:
                jmin=0
                   
            # Elevations
            elevations = elev[imin:imax,jmin:jmax]
            
            if np.sum(elevations) == 0:
                p[i, j], a[i, j],b[i, j], r[i, j] = np.nan, np.nan, np.nan, np.nan
            else:
                           
                # Find 16 neighboring pixels - about ~100 km given that MERRA-2 grid cell size is ~50 x 50 km
                # (40075*np.cos(np.deg2rad(70)))/360
                neighbors = climatology[imin:imax,jmin:jmax]
                
                # Linear regression
                coeffs = stats.linregress(elevations.flatten(), neighbors.flatten())
                
                # Append to grid
                p[i, j], a[i, j],b[i, j], r[i, j] = coeffs[3], coeffs[1], coeffs[0], coeffs[2]
    
    return p, a, b, r

def resample_grids(p, a, b, r, orig):
    
    # Reproject
    p_resample = kd_tree.resample_nearest(source_def, p, target_def, radius_of_influence=50000)
    a_resample = kd_tree.resample_nearest(source_def, a, target_def, radius_of_influence=50000)
    b_resample = kd_tree.resample_nearest(source_def, b, target_def, radius_of_influence=50000)
    r_resample = kd_tree.resample_nearest(source_def, r, target_def, radius_of_influence=50000)
    orig_resample = kd_tree.resample_nearest(source_def, orig, target_def, radius_of_influence=50000)
    
    # Mask
    p_resample[mask == 0] = np.nan
    a_resample[mask == 0] = np.nan
    b_resample[mask == 0] = np.nan
    r_resample[mask == 0] = np.nan
    orig_resample[mask == 0] = np.nan
    
    return p_resample, a_resample, b_resample, r_resample, orig_resample

def downscale(climatology):
    
    # Downscale stats
    p, a, b, r = produce_stats(climatology)
    
    # Resample
    p_resample, a_resample, b_resample, r_resample, orig_resample = resample_grids(p, a, b, r, climatology)

    # Downscale
    downscale = a_resample + (ismip_1km['SRF'].values * b_resample)
    
    # Replace non-significant values to the original grid value
    downscale[p_resample > 0.1] = orig_resample[p_resample > 0.1]

    return downscale, p_resample, r_resample

#%%

swd_allsky_downscaled = np.zeros((2881, 1681, 1))
swd_clrsky_downscaled = np.zeros((2881, 1681, 1))
lwd_allsky_downscaled = np.zeros((2881, 1681, 1))
lwd_clrsky_downscaled = np.zeros((2881, 1681, 1))
t2m_downscaled = np.zeros((2881, 1681, 1))

swd_allsky_p = np.zeros((2881, 1681, 1))
swd_clrsky_p = np.zeros((2881, 1681, 1))
lwd_allsky_p = np.zeros((2881, 1681, 1))
lwd_clrsky_p = np.zeros((2881, 1681, 1))
t2m_p = np.zeros((2881, 1681, 1))


swd_allsky_r = np.zeros((2881, 1681, 1))
swd_clrsky_r = np.zeros((2881, 1681, 1))
lwd_allsky_r = np.zeros((2881, 1681, 1))
lwd_clrsky_r = np.zeros((2881, 1681, 1))
t2m_r = np.zeros((2881, 1681, 1))

for f in range(len(merra_sw_files)):
    print(f)

    # Produce climatology for allsky shortwave
    swd_allsky = np.nanmean(xr.open_dataset(merra_sw_files[f])['swd_allsky'].values, axis=0)
    
    # Downscale daily climatology for allsky shortwave
    swd_allsky_inter, swd_allsky_inter_p, swd_allsky_inter_r = downscale(swd_allsky)
    
    # Produce climatology for allsky shortwave
    swd_clrsky = np.nanmean(xr.open_dataset(merra_sw_files[f])['swd_clrsky'].values, axis=0)
    
    # Downscale daily climatology for allsky shortwave
    swd_clrsky_inter, swd_clrsky_inter_p, swd_clrsky_inter_r = downscale(swd_clrsky)
    
    # Produce climatology for allsky longwave
    lwd_allsky = np.nanmean(xr.open_dataset(merra_lw_files[f])['lwd_allsky'].values, axis=0)
    
    # Downscale climatology for allsky longwave
    lwd_allsky_inter, lwd_allsky_inter_p, lwd_allsky_inter_r = downscale(lwd_allsky)
        
    # Produce climatology for allsky longwave
    lwd_clrsky = np.nanmean(xr.open_dataset(merra_lw_files[f])['lwd_clrsky'].values, axis=0)
    
    # Downscale climatology for allsky longwave
    lwd_clrsky_inter, lwd_clrsky_inter_p, lwd_clrsky_inter_r = downscale(lwd_clrsky)
    
    # Produce climatology for allsky longwave
    t2m = np.nanmean(xr.open_dataset(merra_tm_files[f])['t2m'].values, axis=0)
    
    # Downscale climatology for allsky longwave
    t2m_inter, t2m_inter_p, t2m_inter_r = downscale(t2m)
    
    # Stack
    swd_allsky_downscaled = np.dstack((swd_allsky_downscaled, swd_allsky_inter))
    swd_clrsky_downscaled = np.dstack((swd_clrsky_downscaled, swd_clrsky_inter))
    lwd_allsky_downscaled = np.dstack((lwd_allsky_downscaled, lwd_allsky_inter))
    lwd_clrsky_downscaled = np.dstack((lwd_clrsky_downscaled, lwd_clrsky_inter))
    t2m_downscaled = np.dstack((t2m_downscaled, t2m_inter))
    
    swd_allsky_p = np.dstack((swd_allsky_p, swd_allsky_inter_p))
    swd_clrsky_p = np.dstack((swd_clrsky_p, swd_clrsky_inter_p))
    lwd_allsky_p = np.dstack((lwd_allsky_p, lwd_allsky_inter_p))
    lwd_clrsky_p = np.dstack((lwd_clrsky_p, lwd_clrsky_inter_p))
    t2m_p = np.dstack((t2m_p, t2m_inter_p))
    
    swd_allsky_r = np.dstack((swd_allsky_r, swd_allsky_inter_r))
    swd_clrsky_r = np.dstack((swd_clrsky_r, swd_clrsky_inter_r))
    lwd_allsky_r = np.dstack((lwd_allsky_r, lwd_allsky_inter_r))
    lwd_clrsky_r = np.dstack((lwd_clrsky_r, lwd_clrsky_inter_r))
    t2m_r = np.dstack((t2m_r, t2m_inter_r))


swd_allsky_downscaled = swd_allsky_downscaled[:,:,1:]
swd_clrsky_downscaled = swd_clrsky_downscaled[:,:,1:]
lwd_allsky_downscaled = lwd_allsky_downscaled[:,:,1:]
lwd_clrsky_downscaled = lwd_clrsky_downscaled[:,:,1:]
t2m_downscaled = t2m_downscaled[:,:,1:]

swd_allsky_p = swd_allsky_p[:,:,1:]
swd_clrsky_p = swd_clrsky_p[:,:,1:]
lwd_allsky_p = lwd_allsky_p[:,:,1:]
lwd_clrsky_p = lwd_clrsky_p[:,:,1:]
t2m_p = t2m_p[:,:,1:]

swd_allsky_r = swd_allsky_r[:,:,1:]
swd_clrsky_r = swd_clrsky_r[:,:,1:]
lwd_allsky_r = lwd_allsky_r[:,:,1:]
lwd_clrsky_r = lwd_clrsky_r[:,:,1:]
t2m_r = t2m_r[:,:,1:]

###############################################################################
# Save 1 km dataset to NetCDF
###############################################################################
dataset = netCDF4.Dataset(path + 'allwave-t2m-downscaled.nc', 'w', format='NETCDF4_CLASSIC')
print('Creating %s' % path + 'allwave-t2m-downscaled.nc')
dataset.Title = "Downscaled summer allsky/clearsky downward shortwave/longwave radiation from MERRA-2"
import time
dataset.History = "Created " + time.ctime(time.time())
dataset.Projection = "WGS 84"
dataset.Reference = "Ryan, J. C. et al. (unpublished)"
dataset.Contact = "jryan4@uoregon.edu"
    
# Create new dimensions
lat_dim = dataset.createDimension('y', swd_allsky_downscaled.shape[0])
lon_dim = dataset.createDimension('x', swd_allsky_downscaled.shape[1])
data_dim = dataset.createDimension('z', swd_allsky_downscaled.shape[2])

    
# Define variable types
Y = dataset.createVariable('latitude', np.float32, ('y','x'))
X = dataset.createVariable('longitude', np.float32, ('y','x'))

y = dataset.createVariable('y', np.float32, ('y'))
x = dataset.createVariable('x', np.float32, ('x'))
z = dataset.createVariable('z', np.float32, ('z'))
    
# Define units
Y.units = "degrees"
X.units = "degrees"
   
# Create the actual 3D variable
allsky_swd_nc = dataset.createVariable('swd_allsky', np.float32, ('y','x','z'))
clrsky_swd_nc = dataset.createVariable('swd_clrsky', np.float32, ('y','x','z'))
allsky_lwd_nc = dataset.createVariable('lwd_allsky', np.float32, ('y','x','z'))
clrsky_lwd_nc = dataset.createVariable('lwd_clrsky', np.float32, ('y','x','z'))
t2m_nc = dataset.createVariable('t2m', np.float32, ('y','x','z'))

allsky_swd_p_nc = dataset.createVariable('swd_allsky_p', np.float32, ('y','x','z'))
clrsky_swd_p_nc = dataset.createVariable('swd_clrsky_p', np.float32, ('y','x','z'))
allsky_lwd_p_nc = dataset.createVariable('lwd_allsky_p', np.float32, ('y','x','z'))
clrsky_lwd_p_nc = dataset.createVariable('lwd_clrsky_p', np.float32, ('y','x','z'))
t2m_p_nc = dataset.createVariable('t2m_p', np.float32, ('y','x','z'))

allsky_swd_r_nc = dataset.createVariable('swd_allsky_r', np.float32, ('y','x','z'))
clrsky_swd_r_nc = dataset.createVariable('swd_clrsky_r', np.float32, ('y','x','z'))
allsky_lwd_r_nc = dataset.createVariable('lwd_allsky_r', np.float32, ('y','x','z'))
clrsky_lwd_r_nc = dataset.createVariable('lwd_clrsky_r', np.float32, ('y','x','z'))
t2m_r_nc = dataset.createVariable('t2m_r', np.float32, ('y','x','z'))


# Write data to layers
Y[:] = ismip_1km['lat'].values
X[:] = ismip_1km['lon'].values
y[:] = ismip_1km['lat'].values[:,0]
x[:] = ismip_1km['lon'].values[0,:]
allsky_swd_nc[:] = swd_allsky_downscaled.astype(np.float32)
clrsky_swd_nc[:] = swd_clrsky_downscaled.astype(np.float32)
allsky_lwd_nc[:] = lwd_allsky_downscaled.astype(np.float32)
clrsky_lwd_nc[:] = lwd_clrsky_downscaled.astype(np.float32)
t2m_nc[:] = t2m_downscaled.astype(np.float32)

allsky_swd_p_nc[:] = swd_allsky_p.astype(np.float32)
clrsky_swd_p_nc[:] = swd_clrsky_p.astype(np.float32)
allsky_lwd_p_nc[:] = lwd_allsky_p.astype(np.float32)
clrsky_lwd_p_nc[:] = lwd_clrsky_p.astype(np.float32)
t2m_p_nc[:] = t2m_p.astype(np.float32)

allsky_swd_r_nc[:] = swd_allsky_r.astype(np.float32)
clrsky_swd_r_nc[:] = swd_clrsky_r.astype(np.float32)
allsky_lwd_r_nc[:] = lwd_allsky_r.astype(np.float32)
clrsky_lwd_r_nc[:] = lwd_clrsky_r.astype(np.float32)
t2m_r_nc[:] = t2m_r.astype(np.float32)

z[:] = np.arange(2002,2024)

print('Writing data to %s' % path + 'allwave-t2m-downscaled.nc')
    
# Close dataset
dataset.close()

#%%

