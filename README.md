# mb-methods
Companion GitHub repository to "Bounding the Moment Deficit Rate on Crustal Faults using Geodetic Data: Methods"
by Jeremy Maurer, Paul Segall, and Andrew Bradley. The scripts and functions are set up to replicate the synthetic 
test figures given in the paper. The main files are "mb.m" for the 2-D anti-plane fault and "Coverage_test_3Dfault.m"
for the 3-D simple planar finite faults.   

## 1-D Fault 
mb.m: 
This is a function to reproduce the synthetic tests using the anti-plane strike-slip fault geometry of Fig. 1 in 
the paper (this includes Fig.'s 2, 4, 6a,b, 7, 8, and 11). The file includes function calls that will approximately 
replicate the figures in the paper that use the 1D anti-plane fault geometry. 

## 3-D Fault
Coverage_test_3Dfault.m:
This script calls a series of functions and can replicate Fig.'s 3 and 6c,d. Default values for the parameters are given in
the script, and changes needed to create each plot are listed. 

## License
This code is licensed under the MIT license. It is free to use and distribute. 

