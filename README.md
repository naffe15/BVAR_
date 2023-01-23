# BVAR_
 
Empirical macro toolbox

by F. Ferroni and F. Canova

This repository contains MATLAB functions, routines and documentation to estimate VARs, factor models and local projections with classical or Bayesian methods. The toolbox allows a researcher to conduct inference under various prior assumptions on the parameters, to produce point and density forecasts, to compute spillovers/connectedness across units of a network and to trace out the causal effect of shocks using a number of identification schemes. The toolbox is equipped to handle missing observations, mixed frequencies and time series with large cross-section information (e.g. panels of VAR, factor models and FAVAR). It also contains a number of routines to extract cyclical information and to date business cycles. We describe the methodology employed and implementation of the functions with a number of practical examples.

The matlab 2020 release has a built in function called 'bvar' which causes a crash with previous versions of the toolbox. In current version, we create the function 'bvar_.m' and substitute 'bvar.m' with 'bvar_.m' in all the examples of the tutorial. Codes are backwarad compatible so that for matlab releases earlier than 2020 'bvar.m' still works. 

Link to the youtube tutorial: https://www.youtube.com/channel/UCDepPX4wbdkIqfg438J0h_g 

Matlab toolbox requirements for each example
exercise_1_classical.m: statistics_toolbox; <br />
exercise_2_minn.m:      statistics_toolbox; optimization_toolbox; <br />
exercise_3_irf.m:       statistics_toolbox; optimization_toolbox; <br />
exercise_4_mfvar.m:     statistics_toolbox; optimization_toolbox; <br />
exercise_5_favar.m:     statistics_toolbox; optimization_toolbox; <br />
example_6_VARX.m:     statistics_toolbox; optimization_toolbox; <br />
example_7_LP.m:     statistics_toolbox; optimization_toolbox; <br />
example_8_panels.m:     statistics_toolbox; optimization_toolbox; <br />
example_9_prediction.m:     statistics_toolbox; optimization_toolbox; <br /> 
example_10_VAR_heterosked.m:     statistics_toolbox; optimization_toolbox; <br /> 
example_11_connectedness.m :     statistics_toolbox; optimization_toolbox; <br />
example_12_bdfm.m:     statistics_toolbox; optimization_toolbox; <br />      
