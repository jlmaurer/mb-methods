# mb-methods
Companion GitHub repository to "Bounding the Moment Deficit Rate on Crustal Faults using Geodetic Data: Methods"
Jeremy Maurer, Paul Segall, Andrew Bradley

mb.m: 
This is a function to reproduce the synthetic tests using the anti-plane strike-slip fault geometry of Fig. 1 in the paper. 
The file includes function calls that will approximately replicate the figures in the paper that use the 1D anti-plane fault geometry. 

mcmc.m: 
Very simple Markov Chain Monte Carlo - Metropolis sampler. No adaptive step size, cooling schedule, etc. It is called by mb.m. 
Description of the input and output parameters are in the file. 
