.. _troubleshooting:

Troubleshooting
***************

Compiled here are some common troubleshooting approaches
we would recommend if you are encountering difficulties
with different functionalities of CKBIT.

Model Loading Troubleshooting
-----------------------------
CKBIT is designed to store compiled models so they do 
not need to be recompiled each time they are run. While
this saves significant time, sometimes the model storage
can encounter errors if other operations on the computer
conflict with it. For those cases, we recommend opening
the folder where the script is being run and deleting the
stored model. This will cause the model to be compiled 
again and should solve the problem.

MAP Troubleshooting
-------------------
Due to sensitivity of the initial starting conditions 
and complex models being optimized, the MAP functionality
can encounter challenges. We recommend specifying the
input variable seed with an integer and running the 
function. If this run fails, specify seed with a new
integer. In our experience, the MAP functionality 
usually works within 5 different seed specifications.

MCMC Troubleshooting
--------------------
MCMC runs have convergence checks detailed in the 
Appendix of the publication. When these convergence checks
are showing that the sampled values from the distribution
are not well converged, we recommend increasing the control
criteria like so: control={'adapt_delta':0.99999999, 
'max_treedepth':100}. Please reference Figure 8 of the 
publication to see how to properly implement this control
variable. In our experience, this will solve
most issues experienced with these chemical models.

VI Troubleshooting
------------------
VI runs also have convergence checks detailed in the 
Appendix of the publication. When these convergence
checks are showing that the distributions did not 
converge well, please try specifying more informative
priors and running the model again. We have found that 
for some datasets, the VI estimations will not converge
without more informed prior distributions being implemented.

