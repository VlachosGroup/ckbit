.. _repeatability:

Repeatability
*************

Bayesian inference techniques rely upon stochastic
processes. For example, MCMC sampling collects 
samples from probability distributions in a stochastic
manner that enables us to reconstruct the distributions. 
Each time this sampling occurs, the initialization and
proposals of steps have a stochastic nature that generate
unique sets of samples for each run. However, well
converged runs should produce highly similar sample
statistics (such as mean and standard deviation). Again, 
while these sample statistics are highly similar,
the samples from which they are produced are unique. The
MAP functionality also relies upon stochastic processes,
so these results may slightly vary as well.

Sometimes we want to repeat the exact results between
runs of the software. CKBIT enables this functionality 
through the seed variable for each function. Note, the 
init variable for MAP and MCMC functions must also be set 
for repeat results to be obtained, but this is easy to 
accomplish by not specifying the init variables in the 
input so the default variables are specified. 

Users should note that repeatability occurs on a given
computer, but not across different computers. To clarify,
computer A will always give the same set of results (call
them results A) if the seed and init variables have been
specified in a script. Computer B will always give the same
set of results (call them results B) for that same specified 
script. However, results A will be different from results
B. This is akin to the collected samples of results A being 
different from the collected samples of results B. Again, 
while the samples will be different, the sample statistics
should be highly similar for well converged runs. 
The repeatability is limited to within a given computer 
and does not extend across computers because of the way 
C++ interacts with computer environments. Specifically, C++
is a coding language that does not fully specify how floating 
point arithmetic is implemented. More details on this 
repeatability topic can be found here:
https://mc-stan.org/docs/2_23/reference-manual/reproducibility-chapter.html

While the exact same numerical results cannot be generated
across different computers, the summary results will be 
highly similar for well converged runs. This is the 
standard for a stochastic process like Bayesian inference. 