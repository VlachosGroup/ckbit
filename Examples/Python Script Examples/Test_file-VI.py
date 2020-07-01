# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 12:07:40 2020

@author: Jane
"""

# The VI estimation. This yields estimates of the posterior distributions of   
# each parameter being estimated, but using the VI technique instead of the MCMC.     
# VI is better than MCMC at generating a large number of samples, but is a less      
# robust technique. It is still in its experimental implementation phase.   
# We demonstrate VI estimation below.  

# We can also specify prior distributions and run inference with them. The following   
# example shows the implementation of the follow prior distributions:   
# For the A0 term of rxn 1, a normal distribution with a mean of 11 and standard deviation of 4.   
# For the Ea term of rxn 1, a normal distribution with a mean of 95 and standard deviation of 5.   
# For the Ea term of rxn 1, a normal distribution with a mean of 0.5 and standard deviation of 0.1.   
# For the Ea term of rxn 2, a normal distribution with a mean of 140 and standard deviation of 5.   
# All prior distribution specification must follow Stan's implementation forms:     
# https://mc-stan.org/docs/2_23/functions-reference/unbounded-continuous-distributions.html

#These output values can be compared against the values provided in the attached Word document   
#named "Test File - Outputs Comparison" in the same folder as this script.

from ckbit import cstr

#Import data
file = './CSTR_Data.xlsx'

#Run VI estimation with standard priors
v1, v2 = cstr.VI(filename=file, pH=True, seed=3, output_samples=10000,
                priors = ['A0[1] ~ normal(11,4)',
                          'Ea[1] ~ normal(95,5)',
                          'A0[2] ~ normal(17,4)',
                          'Ea[2] ~ normal(140,5)'])