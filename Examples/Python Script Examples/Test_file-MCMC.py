# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 12:07:40 2020

@author: Jane
"""

# The MCMC estimation. This yields estimates of the posterior distributions of   
# each parameter being estimated.

#These output values can be compared against the values provided in the attached Word document   
#named "Test File - Outputs Comparison" in the same folder as this script.

from ckbit import cstr
import arviz 

#Import data
file = './CSTR_Data.xlsx'

#Run MCMC estimation with standard priors
m1, m2 = cstr.MCMC(filename=file, pH=True, seed=3, iters=10000)

# There are convergence checks to ensure that these samples can be relied upon.   
# These checks are discussed in detail in the published article. This run passes all    
# those checks, and offers a successful inference we can trust.

# It is important to visualize the correlation that exists between the samples of    
# the parameters, which we can accomplish with a pair plot.

#Generate pairplot
arviz.plot_pair(m1)