# Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records.

Tool for generating target spectrum-compatible fully-nonstationary artificial seismic ground motions 

The produced suites of ground motions can match a given target mean spectrum and target variability for the whole period range of interest.

Works also for zero variability, i.e. only mean matching

Further details are provided in the following document:

Yanni H., Fragiadakis M., and Mitseas I.P. Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records. Earthquake Engineering and Structural Dynamics. In review.

Version 1.0 created by Hera Yanni, first release: 28th of April, 2024 

## How to run
* Run the MAIN_generate_acc.m 
*  Download the files and run the main MAIN_generate_acc.m in MATLAB

## Main features:
* The target spectrum can be a smooth code spectrum or a GMM
* Variability can be obtained from a GMM
* Spectral correlation between periods can be taken into account
* There are two variants for the propabilistic generation of the artificial accelerograms

## Additional features included:
* The new time-modulating function proposed in Yanni H., Fragiadakis M., and Mitseas I.P. Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records.
* Time modulating function of Jennings PC, Housner GW, Tsai C. computation 
* Husid function computation
* Heaviside step function computation

## Copyright (c) 2024

Hera Yanni

Structural Engineer NTUA, MSc in ADERS

Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA

email: heragian@mail.ntua.gr, heragian@gmail.com 
