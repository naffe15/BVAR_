%% BVAR tutorial: LOCAL PROJECTIONS
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     27/02/2020, revised  14/10/2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1) classical  LP with  OLS and  IV
% 2) Bayesian LP with calibrated  parameters and  optimal shrinkage
%     both  OLS  and  IV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear;

addpath ../../cmintools/
addpath ../../bvartools/

load DataGK
y = [logip logcpi gs1 ebp];

%% 1) Classical Cholesky
lags               = 12;
options.hor        = 48;
options.conf_sig   = 0.9;
dm1 = directmethods(y,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to gs1)
indx_sho              = 3;   
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
bvarx = bvar_(y,lags,options);
var_irf_sort = sort(bvarx.ir_draws,4);
% add IRF plot
options.add_irfs(:,:,:,1) = var_irf_sort(indx_var,:,indx_sho,round(bvarx.ndraws*0.95));
options.add_irfs(:,:,:,2) = var_irf_sort(indx_var,:,indx_sho,round(bvarx.ndraws*0.05));
% plot LP
plot_irfs_(dm1.ir_lp(indx_var,:,indx_sho,:),options)
pause;

options =rmfield(options,'add_irfs');

%% 2) Classical IV
% load the instruments
[numi,txti,rawi] = xlsread('factor_data.xlsx','factor_data');
% csv is  not read  in  many  matlab stations.
% depends on the country.  transform  csv  into  xls.

% instrument must have the same lenght as the observed data
options.proxy  = nan(length(y),1);
% use the same instrument as GK 
% instruments and data end in 2012m6

% 1: OLS with controls:
% use the instrument/proxy shock directly in the LP
options.proxy(length(T)- length(numi)+1:end) = numi(:,4);
dm2 = directmethods(y,lags,options);

options0 = options;
options0.fontsize = 18;
options0.saveas_strng  = 'IV';
options0.nplots = [2 2];
options0.varnames = {'IP','CPI','1 year rate','EBP'};
% finally, the plotting command
norm = dm2.irproxy_lp(3,1,1,2)*4;
plot_irfs_(dm2.irproxy_lp(:,:,1,:)/norm,options0)

% 2. 2SLS with controls: 
% Regress first the policy varible on the instrument/proxy shock. Use the
% fitted values of the first stage regression (portion of the policy
% variable explained by the instrument/proxy) as the 
FirstStageReg = ols_reg(gs1,[ones(length(gs1),1) options.proxy]);
options2s.proxy  = nan(length(y),1);
options2s.proxy(FirstStageReg.nindex) = FirstStageReg.yfit;
options2s.hor        = 48;
options2s.conf_sig   = 0.9;
dm2s = directmethods(y,lags,options2s);
% finally, the plotting command
norm = dm2s.irproxy_lp(3,1,1,2)*4;
options0.saveas_strng  = 'IV2S';
plot_irfs_(dm2s.irproxy_lp(:,:,1,:)/norm,options0)


% 3. no controls except for lags of the endogenous
options1 = options0;
irf2plot = nan(size(dm2.irproxy_lp));
for vv = 1 :size(y,2)
    dm3 = directmethods (y(:,vv),lags,options);    
    irf2plot(vv,:,1,:) = dm3.irproxy_lp(:,:,1,:)/norm;
end
% dm2 and  dm3  differ  in  the  conditioning  variables (lags of
% all variables,  lags  of  only  own  variables.
options1.saveas_strng  = 'IV_nocontrols';
options1.varnames = options0.varnames;
% finally, the plotting command
plot_irfs_(irf2plot,options1);

pause;


%% 3) Bayesian Local Projection

% run a VAR on presample data
presample = 96; % 8 years of presample
lags      = 12;
BVAR1     = bvar_(y(1:presample,:),lags);
  

% use the VAR estimates to set the priors for the LP
options.priors.name        = 'Conjugate';
% posterior mean of the VAR AR coeff and constant
options.priors.Phi.mean    = mean(BVAR1.Phi_draws,3);
% average variance of the AR coeff and constant
options.priors.Phi.cov     = diag(mean(var(BVAR1.Phi_draws,0,3),2));
% posterior mean of the Covariance of the VAR residuals 
options.priors.Sigma.scale = mean(BVAR1.Sigma_draws,3);
options.priors.Sigma.df    = size(BVAR1.Phi_draws,1)-2;
options.priors.tau         = 0.5*ones(options.hor); % 
options.K                  = 1000;
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

pause;
% close  all

%% 4) Bayesian  LP  with Optimize Shrinkage
 
options.priors.max_tau = 1; % 
options.max_compute    = 1; % fmin search 
bdm_opt                = directmethods(y(presample+1:end,:),lags,options);

% the plotting command
options.saveas_strng  = 'BLPCholeskyOpt';
plot_irfs_(bdm_opt.ir_blp(indx_var,:,indx_sho,:),options)

% the plotting command
options.saveas_strng  = 'BLPIVOpt';
plot_irfs_(bdm_opt.irproxy_blp(indx_var,:,1,:)*0.25,options)

