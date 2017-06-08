# mb-methods
Companion GitHub repository to "Bounding the Moment Deficit Rate on Crustal Faults using Geodetic Data: Methods"
by Jeremy Maurer, Paul Segall, and Andrew Bradley, currently under review in JGR-Solid Earth. The scripts and functions are set up to replicate the synthetic 
test figures given in the paper. The main files are "mb.m" for the 2-D anti-plane fault and "Coverage_test_3Dfault.m"
for the 3-D simple planar finite faults.   

## 1-D Fault 
mb.m: 
This is a function to reproduce the synthetic tests using the anti-plane strike-slip fault geometry of Fig. 1 in 
the paper (this includes Fig.'s 2, 4, 6a,b, 7, 8, and 11). The file includes function calls that will approximately 
replicate the figures in the paper that use the 1D anti-plane fault geometry. mb.m uses the mcmc.m function, also 
provided here, to generate the mcmc pdfs. 

## 3-D Fault
Coverage_test_3Dfault.m:
This script calls a series of functions and can replicate Fig.'s 3 and 6c,d. Default values for the parameters are given in
the script, and changes needed to create each plot are listed. 

## Dependencies
The mb.m function does not require any additional code, but the Coverage_test_3Dfault.m script requires disloc3d.m, a matlab wrapper/implementation of the Okada (1992) displacement solutions for a fault in an elastic half-space. This software is freely available from https://pangea.stanford.edu/cdfm/software.

## License
This code is licensed under the MIT license. It is free to use and distribute. 

## Acknowledgements
The "patchfault.m" function was created by an unknown previous group member of the Crustal Deformation and Fault Mechanics
group at Stanford University. Inspiration for some parts of this code came from other past group members, including Kaj Johnson and Jessica Murray. 

## Examples
The headers of the main files described above include examples for replicating each of the figures in the paper. 
Some examples for the 1-D fault: 

 [out] = mb('simple_bootstrap_test');
 [out] = mb('compare_cobe_coble');

and for the 3-D fault tests, simply run 

"Coverage_test_3Dfault.m"


