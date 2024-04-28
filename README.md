# Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records.

Tool for generating target spectrum-compatible fully-nonstationary artificial seismic ground motions 

The produced suites of ground motions can match a given target mean spectrum and target variability for the whole period range of interest.

Works also for zero variability, i.e. only mean matching

Further details are provided in the following document:

Yanni H., Fragiadakis M., and Mitseas I.P. Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records. Earthquake Engineering and Structural Dynamics. In review.

Version 1.0 created by Hera Yanni, first release: 28th of April, 2024 

## How to run
*  Download the files and run the main MAIN_generate_acc.m in MATLAB
  
* Run the MAIN_generate_acc.m in MATLAB online (no license is needed for MATLAB online)
  by clicking on this link https://matlab.mathworks.com/open/github/v1?repo=HeraYanni/Propabilistic_generation_of_artificial_accelerograms,

  or pressing this button [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=HeraYanni/Propabilistic_generation_of_artificial_accelerograms)

  Needs to login to a MATLAB account, the working interface is this:

  ![image](https://github.com/HeraYanni/Propabilistic_generation_of_artificial_accelerograms/assets/159805439/7dea1e22-4074-4f56-8ba3-9f1e64136188)



## Main features:
* The target spectrum can be a smooth code spectrum or a GMM
* Variability can be obtained from a GMM
* Spectral correlation between periods can be taken into account
* There are two variants for the probabilistic generation of the artificial accelerograms

## Additional features included:
* The new time-modulating function proposed in Yanni H., Fragiadakis M., and Mitseas I.P. Probabilistic generation of hazard-consistent suites of fully non-stationary seismic records.
* Time modulating function of Jennings PC, Housner GW, Tsai C. computation 
* Husid function computation
* Heaviside step function computation


## Copyright Notice

The present software files, created by Hera Yanni, can be freely employed, downloaded, and shared with others as long as proper credit is attributed to the creator:

Hera Yanni (c) 2024

Structural Engineer NTUA, MSc in ADERS

Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA

email: heragian@mail.ntua.gr, heragian@gmail.com 

Please note that there is no warranty or liability of the creator associated with this free software.
