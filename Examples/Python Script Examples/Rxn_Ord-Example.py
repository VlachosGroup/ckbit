# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 14:23:57 2020

@author: Jane
"""

import arviz
from ckbit import rxn_ord

#Import data
file = './RO_data.xlsx'

#Run MAP estimation with standard priors
map1 = rxn_ord.MAP(filename=file)

#Run MCMC estimation with standard priors
m1, m2 = rxn_ord.MCMC(filename=file,control={'adapt_delta':0.99999999, 
                            'max_treedepth':100}, iters=1000, chains=2)

#Generate pairplot
arviz.plot_pair(m1)

#Run VI estimation with standard priors
v1, v2 = rxn_ord.VI(filename=file)

#Process data
data_dict={'intercept':v1['sampler_params'][0], 
           'rxn_ord':v1['sampler_params'][1], 
           'sigma':v1['sampler_params'][2]} 

#Generate pairplot
arviz.plot_pair(data_dict)

#Run MCMC estimation with specified priors
p1, p2 = rxn_ord.MCMC(filename=file,control={'adapt_delta':0.99999999,
                           'max_treedepth':100}, iters=1000,
                            priors = ['rxn_ord ~ normal(1,0.05)'])
