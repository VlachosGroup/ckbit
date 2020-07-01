# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 12:07:40 2020

@author: Jane
"""

# The MAP estimation. This yields point estimates of the modes of the posterior.    
# These estimates are the values that fit the model best given the data and priors.

#These output values can be compared against the values provided in the attached Word document   
#named "Test File - Outputs Comparison" in the same folder as this script.

from ckbit import cstr

#Import data
file = './CSTR_Data.xlsx'

#Run MAP estimation with standard priors
map_vals = cstr.MAP(filename=file, pH=True, seed=3)