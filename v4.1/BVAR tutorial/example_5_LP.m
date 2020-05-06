%% BVAR tutorial: FORECASTS
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc

addpath ../../cmintools/
addpath ../../v4.1/

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
% Compare with BVAR estimates
bvar_ = bvar(y,lags,options);
var_irf_sort = sort(bvar_.ir_draws,4);
% add IRF plot
options.add_irfs(:,:,:,1) = var_irf_sort(indx_var,:,indx_sho,round(bvar_.ndraws*0.95));
options.add_irfs(:,:,:,2) = var_irf_sort(indx_var,:,indx_sho,round(bvar_.ndraws*0.05));
% plot LP
plot_irfs_(dm1.ir_lp(indx_var,:,indx_sho,:),options)

options =rmfield(options,'add_irfs');

%% 2/ IV
% load the instruments
[numi,txti,rawi] = xlsread('factor_data.csv','factor_data');
% instrument must have the same lenght as the observed data
options.proxy  = nan(length(y),1);
% use the same instrument as GK 
% instruments and data end in 2012m6
options.proxy(length(T)- length(numi)+1:end) = numi(:,4);

dm2 = directmethods(y,lags,options);

options0= options;
options0.fontsize =18;
options0.saveas_strng  = 'IV';
options0.nplots = [1 1];
options0.varnames = {'IP','CPI','1 year rate','EBP'};
% finally, the plotting command
norm = dm2.irproxy_lp(3,1,1,2)*4;
plot_irfs_(dm2.irproxy_lp(:,:,1,:)/norm,options0)

options1 = options0;
for vv = 1 :size(y,2)
    dm3 = directmethods (y(:,vv),lags,options);    
    options1.saveas_strng  = ['IV_var' num2str(vv)];
    options1.varnames = options0.varnames(vv);
    % finally, the plotting command
    plot_irfs_(dm3.irproxy_lp(:,:,1,:)/norm,options1)
end

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
options.priors.tau         = 0.5*ones(options.hor); % 

options.proxy(1:presample,:) =[];

bdm = directmethods(y(presample+1:end,:),lags,options);

options.saveas_strng  = 'BLPCholesky';
% the plotting command
options.conf_sig   = 0.9;
options.fontsize   = 12;
plot_irfs_(bdm.ir_blp(indx_var,:,indx_sho,:),options)

options.saveas_strng  = 'BLPIV';

norm = median(bdm.irproxy_blp(3,1,1,:),4)*4;
% the plotting command
plot_irfs_(bdm.irproxy_blp(indx_var,:,1,:)/norm,options)

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


