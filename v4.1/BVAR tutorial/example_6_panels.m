%% BVAR tutorial: IRF with crossection
% Author:   Filippo Ferroni
% Date:     01/05/2020

clear all
close all
clc
addpath ../../cmintools/
addpath ../../v4.1/

%% %=======================================================================
%% POOLING %%%
%% %=======================================================================
%% 0\ Banks (Pooling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating the transmission of a lending  shock  to  lending rates and
% deposit rates  in a VAR assuming bank by bank, then take average (
% consistent estimator when T -> oo

% the data  is a  monthly panel of  100  data  points  for  50  banks.
% load data 
load DataBanks
[T,NBanks]   = size(LendingRate);

lags        = 2;
opt.K       = 1;    % Draws from the posterior not needed
for i  = 1 : NBanks
    % construct the banks i database 
    yi   = [LendingRate(:,i) DepositRate(:,i)];    
    % estimate VAR(2) 
    bvar0 = bvar(yi,lags,opt);
    % store IRFs
    irfs_to_plot(:,:,:,i)  = bvar0.ir_ols;
end

% Customize the IRF plot
% variables names for the plots
options.varnames      = {'Lending Rate','Deposit Rate'};  
% names of the directory where the figure is saved
options.saveas_dir    = './panels_plt';
% names of the figure to save
options.saveas_strng  = 'TSAverageIRF';
% name of the shock
options.shocksnames   = options.varnames;
% additional 90% HPD set
options.conf_sig_2    = 0.95;   
% plot appeareance
options.nplots = [2 2];
% add mean response 
options.add_irfs = mean(irfs_to_plot,4);
% finally, the plotting command
plot_all_irfs_(irfs_to_plot,options);

%% 1\ Banks (Pooling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WG_dynamic.m, fabio canova, version 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating the transmission of a lending  shock  to  lending rates and
% deposit rates  in a VAR assuming  homogenous dynamics and  fixed  effects
% using  a fixed effect  estimator.
% the data  is a  monthly panel of  100  data  points  for  50  banks.

% NOTE: IF THE DATA  HAS DYNAMIC HETEROGENEITY THE  ESTIMATOR COMPUTED HERE
% IS INCONSISTENT
LendingRate_ = reshape(demean(LendingRate),T*NBanks,1);
DepositRate_ = reshape(demean(DepositRate),T*NBanks,1);

ypooled  = [LendingRate_ , DepositRate_];


lags        = 2;
bvar1       = bvar(ypooled,lags);

% IRF to PLOT
irfs_to_plot           = bvar1.ir_draws;
% names of the figure to save
options.saveas_strng  = 'Pooling';
% the plotting command
plot_all_irfs_(irfs_to_plot,options);


%% 2/ Banks (partial pooling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% partial_pool_bayesian.m, fabio canova, version 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimating the transmission of a lending  shock  to  lending  and  deposits
% rates assuming  heterogenous dynamics and  fixed  effects
% the data  is a  monthly panel of  100  data  points  for  50  banks
% Using  Bayesian partial pooling estimator  for  coefficients
% Bpost=(X'*X + (1./gam)*eye(k,k))\(X'*X*Bet+ (1./gam)*barBet);
% barBet =  prior  mean, SigmaV=gam=*Sigma prior dispersion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LendingRateDM = demean(LendingRate);
DepositRateDM = demean(DepositRate);
% number of variables
Nv = 2; 

% Preparing the Pooling Parameters
% prior  cross  sectional  mean
k       = lags*Nv + 1;
% pooling  toward  a  prior  mean  of  zero
barBet  = zeros(k,Nv);  
% pooling  toward  a random  walk
barBet(1,1)         = 1.0;
barBet(lags+1,2)   = 1.0;  
% Pooling Shrinkage
gam = 0.1;          % pooling  parameter :  if  very  large  just  OLS

options.priors.name     = 'Conjugate';
options.priors.Phi.mean = barBet;
options.priors.Phi.cov  = gam * eye(size(barBet,1));
options.conf_sig_2      = 0.9;   

for i  = [1 NBanks]    
    yi  = [LendingRateDM(:,i) DepositRateDM(:,i)];    
    bvar_ = bvar(yi,lags,options);    
    % names of the figure to save
    options.saveas_strng  = ['Bank#' num2str(i)];
    plot_all_irfs_(bvar_.ir_draws,options)

end

%% 3/ Banks (Exchangable prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use a subset of units to form the prior for the remaining units
% Estimating the transmission of a lending  shock  to  lending  and deposits
% rates assuming  heterogenous dynamics and  fixed  effects
% the data  is a  monthly panel of  100  data  points  for  50  banks
% Using  Bayesian partial pooling estimator  for  coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare the Prior
% use the first 20 units to get the pool estimator
N1           = 20; 
LendingRate_ = reshape(demean(LendingRate(:,1:N1)),T*N1,1);
DepositRate_ = reshape(demean(DepositRate(:,1:N1)),T*N1,1);
yp           = [LendingRate_ , DepositRate_];
bvar2        = bvar(yp,lags);


options.priors.name     = 'Conjugate';
options.priors.Phi.mean = mean(bvar2.Phi_draws,3);
% Exchangable prior Shrinkage
% if  very  large  just  OLS, if  very  small perfect  pooling.
gam = 0.1;          
options.priors.Phi.cov  = gam * eye(size(options.priors.Phi.mean,1));
options.K = 1000;
i1        = 0;
for i  = [N1 + 1 : NBanks]    
    i1 = i1 + 1;
    % construct the banks i database 
    yi  = [LendingRate(:,i) DepositRate(:,i)];
    % estimate bivarite VAR
    bvar_e = bvar(yi,lags,options);
    % store IRFs
    irfs_to_plot2(:,:,:,i1)  = bvar_e.ir_ols;
    if i == NBanks
        options.saveas_strng  = ['BayesianPartialPoolingBank#' num2str(i)];
        plot_all_irfs_(bvar_e.ir_draws,options)
    end
end

% names of the figure to save
options.saveas_strng  = 'BayesianPartialPooling';
% the plotting command
plot_all_irfs_(irfs_to_plot2,options);



return
%%=========================================================================
%% Extra\ Countries
% load the data
load DataPooling
% Time:         1978m1 to 2012m8
cnames = {'uk','us','jp','de'};
Nc = length(cnames);
% Varibles:     IPI, CPI, 1Y GOVT YIELD (LTR), Policy Rate (STR)
vnames = {'ipi','cpi','ltr','str'};
Nv = length(vnames);

% Transform the data: take log-differences and demean.
% reshape the data so that 
% [ ipi_uk, cpi_uk, ltr_uk, str_uk;
%   ipi_us, cpi_us, ltr_us, str_us;
%   ipi_jp, cpi_jp, ltr_jp, str_jp;
%   ipi_de, cpi_de, ltr_de, str_de]
T       = size(time,1)-1;
ypooled = nan(T*Nc,Nv);

for v = 1 : Nv % iterate over var
    dta = [];
    for c = 1 : Nc % iterate over countries
        eval(['tmp = ' vnames{v} '_' cnames{v} ';']);
        dta =  [dta ; demean(100*diff(log(tmp)))];
    end
    ypooled(:,v) = dta;
end

lags        = 4 ;
options.hor = 36;
bvar1       = bvar(ypooled,lags,options);

% Define the IRF of Interest
% index of the shocks of interest (shock to str)
indx_sho              = [4];   
% IRF to PLOT
irfs_to_plot           = bvar1.ir_draws(:,:,indx_sho,:);

% Customize the IRF plot
% variables names for the plots
options.varnames      = vnames;  
% names of the directory where the figure is saved
options.saveas_dir    = './panels_plt';
% names of the figure to save
options.saveas_strng  = 'pooling';
% name of the shock
options.shocksnames   = {'MP'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% finally, the plotting command
plot_irfs_(irfs_to_plot,options)

