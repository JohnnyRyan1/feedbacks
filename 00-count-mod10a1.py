#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Count how many MOD10A1 tiles downloaded.

"""

# Import modules
import pandas as pd
import glob
import numpy as np
import os

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'

# Define files
files = sorted(glob.glob(path + 'mod10a1/*'))

# Define years
years = np.arange(2002, 2024, 1)

#%%
start = []
end = []
count = []

for i in years:
    
    modis_files = []
    # Get MODIS tiles
    for f in files:
        
        # Get path and filename seperately 
        infilepath, infilename = os.path.split(f) 
        # Get file name without extension            
        infileshortname, extension = os.path.splitext(infilename)
        
        if infileshortname[9:13]  == str(i):
            modis_files.append(f)
    
    # Get path and filename seperately 
    startpath, startfilename = os.path.split(modis_files[0]) 
    # Get file name without extension            
    startshortname, extension = os.path.splitext(startfilename)
    
    # Get path and filename seperately 
    endpath, endfilename = os.path.split(modis_files[-1]) 
    # Get file name without extension            
    endshortname, extension = os.path.splitext(endfilename)
    
    # Append dates
    start.append(startshortname[13:16])
    end.append(endshortname[13:16])
    count.append(len(modis_files))


#%%

df = pd.DataFrame(list(zip(years, start, end, count)), 
                  columns=['year', 'start', 'end', 'count'])




















