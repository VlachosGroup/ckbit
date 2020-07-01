# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 16:15:20 2020

@author: Jane
"""

import arviz
from ckbit import cstr

#Import data
file = './CSTR_Data.xlsx'

#Run MAP estimation with standard priors
map1 = cstr.MAP(filename=file, pH=True, seed=7)

#Run MCMC estimation with standard priors
m1, m2 = cstr.MCMC(filename=file, pH=True)

#Generate pairplot
arviz.plot_pair(m1)

#Run VI estimation with standard priors
v1, v2 = cstr.VI(filename=file, pH=True,
                priors = ['A0[1] ~ normal(11,4)',
                          'Ea[1] ~ normal(95,5)',
                          'A0[2] ~ normal(17,4)',
                          'Ea[2] ~ normal(140,5)'])
