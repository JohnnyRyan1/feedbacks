#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" 

DESCRIPTION

1. Interpolate missing data gaps in MODIS albedo data.

"""

# Import modules
import xarray as xr
import netCDF4
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os

#%%

# Define user
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Import ISMIP 1 km grid
ismip_1km = xr.open_dataset(path + '1km-ISMIP6-GIMP.nc')

# Import MODIS files
modis_files = sorted(glob.glob(path + 'modis-albedo/*.nc'))

# Define destination
dest = path + 'modis-albedo-filled/'

#%%

# Define mask
mask = ismip_1km['GIMP'].values

# Define 3D mask
mask3d = np.repeat(mask[:,:,np.newaxis], 125, axis=2)

#%%

"""
Code for testing single pixels

"""

# Import MODIS data
modis = xr.open_dataset(modis_files[0])

# Find random pixels contained in ice mask
pixels = np.array((np.where(ismip_1km['GIMP'])[0], np.where(ismip_1km['GIMP'])[1]))
random = pixels[:,::5000]

valid = []

for i in range(random.shape[1]):
    # Find pixel that is in ice mask
    point = modis['albedo'][random[0,i], random[1,i], :]
    
    # Replace zeros with nan
    point = point.astype(np.float32)
    point = point.where(point != 0)
    
    # Perform filter (i.e. remove if two standard deviations from mean)
    rolling_mean = point.rolling(z=11, min_periods=1, center=True).mean()
    rolling_std = point.rolling(z=11, min_periods=1, center=True).std()
    rolling_std = rolling_std * 2
    
    # Calculate difference between pixel value and rolling mean
    difference = np.abs(point - rolling_mean)
    
    # Mask values that are more than two standard deviations from the mean
    mask = (difference < rolling_std)
    
    # Calculate 11-day rolling median to be used as the timeseries
    rolling_median = point.where(mask == True).rolling(z=11, min_periods=3, center=True).median()
    
    # Linearly interpolate between values
    linear_interp = rolling_median.interpolate_na(dim="z", method="linear")

    # Optional plot export
    plt.scatter(modis['z'], point)
    plt.scatter(modis['z'], linear_interp)
    plt.axvline(x=18, ls='dashed', color='k')
    plt.axvline(x=110, ls='dashed', color='k')
    plt.savefig(path + 'filling-figures/plot_' + str(random[0,i]) + 
                '_' + str(random[1,i]) + '.png')
    plt.close()
    

#%%

"""

Problem 1:
    Some pixels e.g. (2208, 1247) do not have enought points so linear interpolation
    only provides values for about 15-20 days.

Problem 2:
    Some pixels e.g. (2410, 1186) there are no values at all.

Problem 3:
    Some pixels e.g. (381, 861) there are not enough at the start.

Solution: 
    Since there are are actually quite a lot of pixels without 92 values (~40%),
    the first step should be to replace the NaNs with values from neighboring 
    grid cells. Hopefully most NaNs will be replaced using this method. 
    
    Note that we must ensure that (1) we compute SWnet for all the ice sheet 
    and (2) we have a consistent number of pixels for every year.


"""

#%%

def save2netcdf(dest, file, lats, lons, filled_albedo, year):
    
    ###############################################################################
    # Save 1 km dataset to NetCDF
    ###############################################################################
    
    # Get path and filename seperately 
    infilepath, infilename = os.path.split(file) 
    # Get file name without extension            
    infileshortname, extension = os.path.splitext(infilename)
        
    dataset = netCDF4.Dataset(dest + infilename, 'w', format='NETCDF4_CLASSIC')
    print('Creating %s' % dest + infilename)
    dataset.Title = "Daily gap-filled albedo for %s from MOD10A1 product" %(str(year))
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C. et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', filled_albedo.shape[0])
    lon_dim = dataset.createDimension('x', filled_albedo.shape[1])
    data_dim = dataset.createDimension('z', filled_albedo.shape[2])

        
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
    filled_albedo_nc = dataset.createVariable('filled_albedo', np.int8, ('y','x','z'))

    # Write data to layers
    Y[:] = lats
    X[:] = lons
    x[:] = lons[0,:]
    y[:] = lats[:,0]
    filled_albedo_nc[:] = filled_albedo.astype(np.int8)
    z[:] = np.arange(1,93)
    
    print('Writing data to %s' % dest + infilename)
        
    # Close dataset
    dataset.close()


#%%

for f in modis_files:
    
    # Get path and filename seperately 
    infilepath, infilename = os.path.split(f) 
    # Get file name without extension            
    infileshortname, extension = os.path.splitext(infilename)
    print('Processing... %s' % f)

# =============================================================================
#     if os.path.exists(dest + 'mod10a1-albedo-' + infileshortname[-4:] + '.nc'):
#         print('Skipping... %s' %(dest + 'mod10a1-albedo-' + infileshortname[-4:] + '.nc'))
#     else:
#     
# =============================================================================
    # Import MODIS data
    mod = xr.open_dataset(f)
    
    # Define albedo
    albedo = mod['albedo']
    
    # Mask out ice sheet
    albedo = mod['albedo'].where(mask3d == True)
    
    # Split into 100 x 100
    chunks = np.arange(0, 3000, 100)

    # Create empty arrays for new data
    new_data = np.zeros(albedo.shape).astype(np.int8)
      
    counts = []
    
    for i in range(len(chunks) - 1):
        for j in range(len(chunks) - 1):           
            # Define chunk mask
            valid = albedo[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]]
            
            if np.sum(valid) > 0:
              
                array = albedo[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]]

                # Convert to float
                array = array.astype(np.float32)
                
                # Replace zeros with NaNs
                array = array.where(array != 0)
                
                # Perform filter (i.e. remove if two standard deviations from mean)
                rolling_mean = array.rolling(z=11, min_periods=1, center=True).mean()
                rolling_std = array.rolling(z=11, min_periods=1, center=True).std()
                rolling_std = rolling_std * 2
                
                # Calculate difference between pixel value and rolling mean
                difference = np.abs(array - rolling_mean)
                
                # Mask values that are more than two standard deviations from the mean
                mask = (difference < rolling_std)
                
                # Calculate 11-day rolling median to be used as the timeseries
                rolling_median = array.where(mask == True).rolling(z=11, min_periods=3, center=True).median()
                
                # Linearly interpolate between values
                linear_interp = rolling_median.interpolate_na(dim="z", method="linear")
                
                # Add to new data array
                new_data[chunks[i]:chunks[i+1],chunks[j]:chunks[j+1]] = linear_interp.values
                
                # Compute number of filled values
                counts.append((np.sum(np.isfinite(array.values)),
                               np.sum(np.isfinite(array.values)) - np.sum(np.isfinite(array.where(mask == True).values)),
                               np.sum(np.isfinite(linear_interp.values))))
                 
# =============================================================================
#     # Save as NetCDF
#     save2netcdf(dest, f, mod['latitude'].values, mod['longitude'].values, 
#                 new_data[:,:,18:110], infileshortname[-4:])
# =============================================================================
    
    # Save DataFrame
    df = pd.DataFrame(counts, columns=['original', 'cut', 'interpolated'])
    df.to_csv(path + 'filling-stats/' + infileshortname + '.csv')
        
    
#%%
    
    
    
    
    
    
    
    
    
