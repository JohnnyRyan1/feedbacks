#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

DESCRIPTION

1. Read and stack MOD10A1 HDFs

"""

# Import modules
import numpy as np
import netCDF4
from pyhdf.SD import SD, SDC
import os
import glob
from pyproj import Proj

#%%

# Define user 
user = 'jryan4'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/'

# Define destination to save
dest = path + 'data/modis-albedo-intermediate/'

# Define location of MODIS data
modis_files = sorted(glob.glob(path + 'data/mod10a1/*.hdf'))
print('Number of files in folder = %.0f' %len(modis_files))

# Get number of unique tiles
modis_tiles = []
for file in modis_files:
    modis_tiles.append(os.path.basename(file)[18:23])

tiles = np.unique(np.array(modis_tiles))
years = np.arange(2022, 2024, 1)

tiles=['15v01','15v02','16v00', '16v01', '16v02', '17v00', '17v01', '17v02']

#%%

def save2netcdf(lats, lons, doy, data1, data2, data3, tile, year):
    dataset = netCDF4.Dataset(dest + 'MOD10A1_Albedo_' + tile + '_' + year + '.nc', 
                              'w', format='NETCDF4_CLASSIC')
    #print('Creating %s' %dest + 'MOD10A1_Albedo_' + tile + '_' + year + '.nc')
    dataset.Title = "Snow albedo for tile %s for %s from MOD10A1 product" %(tile, year)
    import time
    dataset.History = "Created " + time.ctime(time.time())
    dataset.Projection = "WGS 84"
    dataset.Reference = "Ryan, J. C., et al. (unpublished)"
    dataset.Contact = "jryan4@uoregon.edu"
        
    # Create new dimensions
    lat_dim = dataset.createDimension('y', data1.shape[0])
    lon_dim = dataset.createDimension('x', data1.shape[1])
    data_dim = dataset.createDimension('z', data1.shape[2])
        
    # Define variable types
    Y = dataset.createVariable('latitude', np.float64, ('y','x'))
    X = dataset.createVariable('longitude', np.float64, ('y','x'))
    
    y = dataset.createVariable('y', np.float64, ('y'))
    x = dataset.createVariable('x', np.float64, ('x'))
    z = dataset.createVariable('z', np.float64, ('z'))
        
    # Define units
    Y.units = "degrees"
    X.units = "degrees"
       
    # Create the actual 3D variable
    albedo = dataset.createVariable('albedo', np.int8, ('y','x','z'))
    qa = dataset.createVariable('quality_flag', np.int8, ('y','x','z'))
    mask = dataset.createVariable('mask', np.int32, ('y','x'))
    day_of_year = dataset.createVariable('day_of_year', np.int16, ('z'))
    
    albedo.units = "unitless"
    
    # Write data to layers
    Y[:] = lats
    X[:] = lons
    x[:] = lons[0,:]
    y[:] = lats[:,0]
    z[:] = doy
    day_of_year[:] = doy
    albedo[:] = data1
    qa[:] = data2
    mask[:] = data3

    #print('Writing data to %s' %dest + 'MOD10A1_Albedo_' + tile + '_' + year + '.nc')
        
    # Close dataset
    dataset.close()
    
#%%

# Loop through tiles
for i in tiles:
    # Get MODIS files
    modis_list = []
    for file in modis_files:
        if os.path.basename(file)[18:23]  == i:
            modis_list.append(file)
    
    print('Processing tile... %s' %i)
    
    # Define lat and lons of tile
    grid_cell = 463.31271653
    upper_left_grid = (-20015109.354, 10007554.677)
    lower_right_grid = (20015109.354, -10007554.677)
    
    upper_left_corner = (upper_left_grid[0] + (2400*int(i[0:2])*grid_cell),
                         upper_left_grid[1] - (2400*int(i[3:])*grid_cell))
    lower_right_corner = (upper_left_corner[0] + (2400*grid_cell),
                         (upper_left_corner[1] - (2400*grid_cell)))
    
    x = np.linspace(upper_left_corner[0], lower_right_corner[0], 2400)
    y = np.linspace(upper_left_corner[1], lower_right_corner[1], 2400)
    
    # Produce grid
    xv, yv = np.meshgrid(x, y)
    
    # Define MODIS grid in pyproj
    modis_grid = Proj('+proj=sinu +R=6371007.181 +nadgrids=@null +wktext')
    
    # Convert to lat, lons
    lon, lat = modis_grid(xv, yv, inverse=True)
    
    for p in years:
        # Get MODIS files
        modis_list_by_years = []
        for file in modis_list:
            if (os.path.basename(file)[9:13] == str(p)) & (int(os.path.basename(file)[13:16]) > 151) &  (int(os.path.basename(file)[13:16]) < 244):
                modis_list_by_years.append(file)
        
        #print('Processing year... %s' %str(p))
        if os.path.exists(dest + 'MOD10A1_Albedo_' + i + '_' + str(p) + '.nc'):
            print('Skipping... %s' %(dest + 'MOD10A1_Albedo_' + i + '_' + str(p) + '.nc'))
        else:
            # Define empty arrays for data
            snow_albedo_grid = []
            quality_grid = []
            doy = []
                 
            for h in range(len(modis_list_by_years)):
                # Get path and filename seperately 
                infilepath, infilename = os.path.split(modis_list_by_years[h]) 
                # Get file name without extension            
                infileshortname, extension = os.path.splitext(infilename)
                
                # Append day of year to list
                doy.append(int(infileshortname[13:16]))
            
                # Read MODIS file
                f = SD(modis_list_by_years[h], SDC.READ)
                
                # Get datasets
                sds_snow_albedo = f.select('Snow_Albedo_Daily_Tile')
                snow_albedo = sds_snow_albedo.get()
                sds_basic_qa = f.select('NDSI_Snow_Cover_Basic_QA')
                basic_qa = sds_basic_qa.get()
                
                # Stack to empty array
                snow_albedo_grid.append(snow_albedo)
                quality_grid.append(basic_qa)
                
                #print('Processing... %.0f out of %.0f' %(h+1, len(modis_list_by_years)))
                       
                # Optional list datasets #
                #datasets_dic = f.datasets()
            
                #for idx,sds in enumerate(datasets_dic.keys()):
                #    print idx,sds
            
                # Optional export csv to test
                #data = d.variables['albedo'][:,:,0]
                #mask = np.ravel(data < 100)
                #df = pd.DataFrame(list(zip(np.ravel(d.variables['x'][:])[mask], 
                #                           np.ravel(d.variables['y'][:])[mask], np.ravel(data)[mask])))
                #df.columns = ['X', 'Y', 'Z']
                #df.to_csv('/media/johnny/RyanLab/jryan11/PolarSnow/Data/MODIS/test.csv')
            
            # Count good values for each pixel for land mask
            data = np.dstack(snow_albedo_grid)
            data[data == 125] = 1
            data[data > 100] = 0
            data[data < 0] = 0
            data_sum = np.sum(data, axis=2)
    
            # Concatenate array and save to netcdf
            save2netcdf(lat, lon, doy, np.dstack(snow_albedo_grid), 
                        np.dstack(quality_grid), data_sum, i, str(p))
            
            
            