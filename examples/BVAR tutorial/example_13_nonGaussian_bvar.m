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
% plotting the empirical distribution of the residuals against the normal
% distribution (--- line). Departures from the dotted lines indicate
% departures from normality. The figure also reports the estimated skewness
% and kurtosis.
opts.varnames      = {'IP','CPI','1 year rate','EBP'};  
plot_statistics_all(bvar0.white_e_ols,opts)


%% run the BVAR with fat-tails (kurtosis)
options.hor = 24;
% options.robust_bayes = 1; % only kurtosis
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
options.shocksnames   = {'EB'};  
% additional 90% HPD set
options.conf_sig_2    = 0.9;   
% add the 90 HPD sets with Normal errors
var_irf_sort = sort(bvar1.ir_draws,4);
options.add_irfs(:,:,:,1) = var_irf_sort(indx_var,:,indx_sho,round(bvar1.ndraws*0.95));
options.add_irfs(:,:,:,2) = var_irf_sort(indx_var,:,indx_sho,round(bvar1.ndraws*0.05));
options.normz = 0;

% finally, the plotting command
plot_irfs_(irfs_to_plot,options)

%%
% For the same identification scheme (choleski), differences across the IRFs
% distributions are small. The posterior distribution of the AR matrices
% does not depend on the distribution of the error term. Only the posterior
% distribution of Sigma does.      

for hh = 1 : bvar1.ndraws
    sig1(:,hh) = vech(bvar1.Sigma_lower_chol_draw(:,:,hh));
    sig0(:,hh) = vech(bvar0.Sigma_lower_chol_draw(:,:,hh));
end

figure('Name','Posterior Distribution of the elements of chol(Sigma)')
for ii = 1 : size(sig0,1)
    subplot(3,4,ii)
    histogram(sig0(ii,:));
    hold on
    histogram(sig1(ii,:));
end    

%%
