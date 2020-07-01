# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 13:18:21 2020

@author: Jane
"""

import arviz
from ckbit import app_ea

#Import data
file = './App_Ea_Data.xlsx'

#Run MAP estimation with standard priors
map1 = app_ea.MAP(filename=file)

#Run MCMC estimation with standard priors
m1, m2 = app_ea.MCMC(filename=file,control={'adapt_delta':0.99999999, 'max_treedepth':100}, 
                      iters=1000, chains=2)

#Generate pairplot
arviz.plot_pair(m1)

#Run VI estimation with standard priors
v1, v2 = app_ea.VI(filename=file)

#Process data
data_dict={'intercept':v1['sampler_params'][0], 
           'Ea':v1['sampler_params'][1], 
           'sigma':v1['sampler_params'][2]} 

#Generate pairplot
arviz.plot_pair(data_dict)

#Run MCMC estimation with specified priors
p1, p2 = app_ea.MCMC(filename=file,control={'adapt_delta':0.99999999,
                           'max_treedepth':100}, iters=1000,
                            priors = ['app_ea ~ normal(90,5)'])