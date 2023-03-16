# BVAR_
 
Empirical macro toolbox

by F. Ferroni and F. Canova

This repository contains MATLAB functions, routines and documentation to estimate VARs, factor models and local projections with classical or Bayesian methods. The toolbox allows a researcher to conduct inference under various prior assumptions on the parameters, to produce point and density forecasts, to compute spillovers/connectedness across units of a network and to trace out the causal effect of shocks using a number of identification schemes. The toolbox is equipped to handle missing observations, mixed frequencies and time series with large cross-section information (e.g. panels of VAR, factor models and FAVAR). It also contains a number of routines to extract cyclical information and to date business cycles. We describe the methodology employed and implementation of the functions with a number of practical examples.

The matlab 2020 release has a built in function called 'bvar' which causes a crash with previous versions of the toolbox. In current version, we create the function 'bvar_.m' and substitute 'bvar.m' with 'bvar_.m' in all the examples of the tutorial. Codes are backwarad compatible so that for matlab releases earlier than 2020 'bvar.m' still works. 

With Matlab 2022b version, the toolbox crashes when it saves figures in 'eps' or 'pdf'. Until the problem is resolved, the tooolbox does not save figures in 'pdf' or 'eps' if uses the Matlab 2022b version.

Link to the youtube tutorial: https://www.youtube.com/channel/UCDepPX4wbdkIqfg438J0h_g 

Matlab toolbox requirements BVAR tutorial: statistics_toolbox; optimization_toolbox; <br />
Matlab toolbox requirements Trend-Cycle-Dating tutorial tutorial: statistics_toolbox; optimization_toolbox; signal_toolbox; <br />
