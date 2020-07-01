# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 15:42:04 2020

@author: Jane
"""

import arviz
from ckbit import pfr

#Import data
file = './PFR_Data.xlsx'

#Run MAP estimation with standard priors
map1 = pfr.MAP(filename=file, pH=True, rel_tol=5E-6, abs_tol=5E-6, max_num_steps=1000)

#Run MCMC estimation with standard priors
m1, m2 = pfr.MCMC(filename=file, pH=True)

#Generate pairplot
arviz.plot_pair(m1)

#Run MCMC estimation with specified priors
p1, p2 = pfr.MCMC(filename=file, pH=True,
                            priors = ['A0[1] ~ normal(10,5)',
                                      'Ea[1] ~ normal(100,5)'])