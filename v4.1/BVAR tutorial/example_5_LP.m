%% BVAR tutorial: FORECASTS
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc

addpath ..\..\cmintools\
addpath ..\..\v4.1\

%% %=========================================================================
%%% DIRECT METHODS %%%
%%=========================================================================

load DataGK
y = [logip logcpi gs1 ebp];

%% 1/ Cholesky

lags               = 12;
options.hor        = 48;
options.conf_sig   = 0.9;
dm1 = directmethods(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = [3];   
% Data order
% 1. logip; 2. logcpi; 3. gs1; 4. ebp
% Change the order of the variables for the plot 
% 1. gs1; 2. logcpi; 3. logip; 4. ebp
indx_var              = [3, 2, 1, 4];

% PLOT IRF 
% Customize the IRF plot
% variables names for the plots
options.varnames      = {'1 year rate','CPI','IP','EBP'};  
% names of the directory where the figure is saved
options.saveas_dir    = './dm_plt';
% names of the figure to save
options.saveas_strng  = 'Cholesky';
% name of the shock
options.shocksnames   = {'MP'};  
% finally, the plotting command
plot_irfs_(dm1.ir_lp(indx_var,:,indx_sho,:),options)

%% 2/ IV

% load the instruments
[numi,txti,rawi] = xlsread('factor_data.csv','factor_data');
% instrument must have the same lenght as the observed data
options.proxy  = nan(length(y),1);
% use the same instrument as GK 
% instruments and data end in 2012m6
options.proxy(length(T)- length(numi)+1:end) = numi(:,4);

dm2 = directmethods(y,lags,options);

options.saveas_strng  = 'IV';
% finally, the plotting command
plot_irfs_(dm2.irproxy_lp(indx_var,:,1,:)*0.25,options)


%% 3/ Bayesian LP 

% run a VAR on presample data
presample = 96; % 8 years of presample
lags      = 12;
bvar_     = bvar(y(1:presample,:),lags);

% use the VAR estimates to set the priors for the LP
options.priors.name        = 'Conjugate';
% posterior mean of the VAR AR coeff and constant
options.priors.Phi.mean    = mean(bvar_.Phi_draws,3);
% average variance of the AR coeff and constant
options.priors.Phi.cov     = diag(mean(var(bvar_.Phi_draws,0,3),2));
% posterior mean of the Covariance of the VAR residuals 
options.priors.Sigma.scale = mean(bvar_.Sigma_draws,3);
options.priors.Sigma.df    = size(bvar_.Phi_draws,1)-2;

options.proxy(1:presample,:) =[];

bdm = directmethods(y(presample+1:end,:),lags,options);

options.saveas_strng  = 'BLPCholesky';
% the plotting command
options.conf_sig   = 0.68;
options.conf_sig_2 = 0.9;
plot_irfs_(bdm.ir_blp(indx_var,:,indx_sho,:),options)

options.saveas_strng  = 'BLPIV';
% the plotting command
plot_irfs_(bdm.irproxy_blp(indx_var,:,1,:)*0.25,options)

%% 4/ Optimize Shrinkage

options.priors.max_tau = 1; % 
options.max_compute    = 1; % fmin search 
bdm_opt                = directmethods(y(presample+1:end,:),lags,options);

% the plotting command
options.saveas_strng  = 'BLPCholeskyOpt';
plot_irfs_(bdm_opt.ir_blp(indx_var,:,indx_sho,:),options)

% the plotting command
options.saveas_strng  = 'BLPIVOpt';
plot_irfs_(bdm_opt.irproxy_blp(indx_var,:,1,:)*0.25,options)


