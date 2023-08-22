%% BVAR tutorial: Non-Gaussian VAR and Indentification of non-Gaussian shocks
% Author:   Filippo Ferroni and Fabio Canova
% Date:     8/6/2023
% Identification based on the paper "Identification Using Higher-Order
% Moments Restrictions" (2023) by P. Andrade, Ferroni F. and L. Melosi 
% https://www.chicagofed.org/publications/working-papers/2023/2023-28


close all; clc; clear;
addpath ../../cmintools/
addpath ../../bvartools/

% load the data
load DataGK
y = [logip logcpi gs1 ebp];

% run the BVAR
lags        = 12;
bvar0       = bvar_(y,lags);

%% plot statistics of the OLS orthogonalized residuals

opts.varnames      = {'IP','CPI','1 year rate','EBP'};  
plot_statistics_all(bvar0.white_e_ols,opts)


%% run the BVAR with fat-tails (kurtosis)
options.hor = 24;
options.robust_bayes = 1; % only kurtosis
options.robust_bayes = 2; % skewness and kurtosis
bvar1       = bvar_(y,lags,options);

%% 1) Cholesky identification of spread shock

% index of the shocks of interest (shock to EBP)
indx_sho              = 4;   
indx_var              = [1 : 4];
% IRF to PLOT
irfs_to_plot           = bvar0.ir_draws(indx_var,:,indx_sho,:);

% Customize the IRF plot
% variables names for the plots
options.varnames      = {'IP','CPI','1 year rate','EBP'};  
% names of the directory where the figure is saved
options.saveas_dir    = './nonguassian_bvar_plt';
% names of the figure to save
options.saveas_strng  = 'cholesky-robust';
% name of the shock
options.shocksnames   = {'MP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)

%% 2) identification using higher-order moments

% TBA