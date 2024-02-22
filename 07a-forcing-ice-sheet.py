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
pos = pd.read_csv(path + 'positive_forcing_results.csv')
neg = pd.read_csv(path + 'negatiave_forcing_results.csv')
clim = pd.read_csv(path + 'ice-sheet-climate.csv')


#%%

# Melt-albedo feedbacks increase SWnet by...
total_albedo = (pos['observed_albedo'] - pos['fixed_snow_all']) / pos['fixed_snow_all']
print(total_albedo.mean())


print(pos['snowline_forcing'].mean())
print(pos['snow_forcing'].mean())
print(pos['ice_forcing'].mean())


