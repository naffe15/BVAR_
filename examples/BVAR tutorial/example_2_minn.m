% example_2_min_new.m VAR inference with Minnesota Prior
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     27/02/2020, revised  14/10/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute  posterior estimate with 
% 1) mixed calibrated/estimated Minnesota prior 
% 2) optimally chosen  Minnesota prior.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
warning off; clear; close all; clc;

addpath ../../cmintools/
addpath ../../bvartools/


% load Euro data
load Data
y= [IPI HICP CORE Euribor1Y M3 EXRATE];% the variables for the VAR

%% Case 1: Estimation with calibrated Minnesota Prior

lags = 6; 
options.max_minn_hyper  = 1;       % start the  optimization routine
BVAR0                    = bvar_(y,lags,options);
pause;


%% Case 2: Estimation with calibrated/ estimated Minnesota Prior
lags = 6; 
options.max_minn_hyper  = 1;       % start the  optimization routine
options.minn_prior_tau  = 10;      % set tau 
options.index_est       = [3 4];   % hyper-parameters to maximize
options.lb              = [0 0];   % the lower bounds 
options.ub              = [20 20]; % the upper bounds
options.max_compute     = 7;       % optimization  by Matlab Simplex
BVAR1                   = bvar_(y,lags,options);
pause;

%% Case 3:  Optimization of one dimension of the Minnesota Prior
%           (without  posterior draws to speed up computations) 
%%  3.1 maximization of  tau
clear options
lags = 6;
% setting the default values for the hyperparameters
hyperpara(1)    = 3;		  % tau
hyperpara(2)    = 0.5;		  % decay
hyperpara(3)    = 5;		  % lambda
hyperpara(4)    = 2;		  % mu
hyperpara(5)    = 2;		  % omega
% setting the options
options.index_est	   = 1:1;    % hyper-parameter over which maximize
options.max_compute    = 2;      % maximize  using Matlab fmincon function
options.lb             = 0.8;    % Lower bound
options.ub             = 10;     % Upper bound
[postmode,logmlike,~] = bvar_max_hyper(hyperpara,y,lags,options);
pause;

%% 3.2    Take  optimal value for hyperparameter(1) and compute optimal 
%%         values for tau, decay, lambda (without  posterior  draws) 
hyperpara(1)            = postmode(1); % use as starting value previous mode
options.index_est       = 1:3;         % set hyper-parameters over which maximize
options.lb              = [0.1 0.1 0.1]; % Lower bounds
options.ub              = [10 10 10];    % Upper bounds
[postmode1,log_dnsty,~] = bvar_max_hyper(hyperpara,y,lags,options);

%% 3.3  Take  optimal value for hyperparameter(1:3) and compute optimal 
%%           values for tau, decay, lambda, mu (without  posterior  draws)

hyperpara(1:3)          = postmode1(1:3); % use as starting value previous mode
options.index_est       = 1:4;         % set hyper-parameters over which maximize
options.lb              = [0.1 0.1 0.1 0.1]; % Lower bounds
options.ub              = [10 10 10 10];    % Upper bounds
[postmode,log_dnsty1,~] = bvar_max_hyper(hyperpara,y,lags,options);
% then run BVAR with  optimal  parameters
BVAR2                   = bvar_(y,lags,options);

%Plotting cholesky responses  to  Monetary  policy  shocks
% Define the IRF of Interest
indx_sho              = 4;   
% Change the order of the variables for the plot 
indx_var              = [4, 1, 2, 3];

% % IRFs to PLOT: compare different max hyperparam
mltple_irfs_to_plot_all(:,:,1,:) = BVAR0.ir_draws(indx_var,:,indx_sho,:);
mltple_irfs_to_plot_all(:,:,2,:) = BVAR1.ir_draws(indx_var,:,indx_sho,:);
mltple_irfs_to_plot_all(:,:,3,:) = BVAR2.ir_draws(indx_var,:,indx_sho,:);

% Customize the plot
% variables names for the plots
options.varnames      = {'1 year rate','IP','HICP','CORE'};  
% name of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% name of the figure to save
options.saveas_strng  = 'diff_minn_opt';
% name of the shock
options.shocksnames   = {'MP opt1','MP opt2','MP opt3'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% the plotting command
plot_all_irfs_(mltple_irfs_to_plot_all,options)

