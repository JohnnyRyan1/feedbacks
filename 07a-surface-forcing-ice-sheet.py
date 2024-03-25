#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Analysis at the ice sheet scale

"""

# Import modules
import pandas as pd


#%%

# Define user
user = 'johnnyryan'

# Define path
path = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/data/'
savepath = '/Users/' + user + '/Dropbox (University of Oregon)/research/feedbacks/figures/'

#%%

# Imoport data
pos = pd.read_csv(path + 'positive-forcing-results.csv')

#%%

# Mean observed SWnet for the entire ice sheet between 2002-2023 
print(pos['observed_albedo'].mean(), '+/-', pos['observed_albedo'].std())


# Mean SWnet for ice sheet with 0.83 albedo between 2002-2023 
print(pos['fixed_snow_all'].mean(), '+/-', pos['fixed_snow_all'].std())

# Contribution of surface albedo change to SWnet
print(((pos['observed_albedo'] - pos['fixed_snow_all']) / pos['fixed_snow_all']).mean())
print(((pos['observed_albedo'] - pos['fixed_snow_all']) / pos['fixed_snow_all']).std())

# Radiative forcing due to snowline fluctuations
print(pos['snowline_forcing'].mean(), '+/-', pos['snowline_forcing'].std())

# Radiative forcing due to snow albedo change
print(pos['snow_forcing'].mean(), '+/-', pos['snow_forcing'].std())

# Radiative forcing due to ice albedo change
print(pos['ice_forcing'].mean(), '+/-', pos['ice_forcing'].std())

# Contribution of surface albedo to SWnet
print(pos['observed_albedo'].mean() - pos['fixed_snow_all'].mean())

#%%










