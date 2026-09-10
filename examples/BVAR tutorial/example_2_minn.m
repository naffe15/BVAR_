% example_2_minn.m VAR inference with Minnesota Prior
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     27/02/2020, revised  19/02/2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute  posterior estimates with 
% 1) mixed calibrated/estimated Minnesota prior hyperparameters
% 2) optimally chosen  Minnesota prior hyperparameters
% 3) Compare  IRFS to  MP shock with  three  hyperparameter choices
% 4) Comparison of log marginal data density and IRFs to an EBP shock 
% between Pandemic Priors and the classic Minnesota Prior
% 5) Minnesota with the long run prior and comparison IRFS to a GDP shock
% with and without it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off; clear; close all; clc;

addpath ../../cmintools/
addpath ../../bvartools/


% load Euro data
load Data
y= [IPI HICP CORE Euribor1Y M3 EXRATE];% the variables for the VAR

%% Case 1: Estimation with cherry-picked Minnesota Prior
lags = 6; 
options.max_minn_hyper  = 1;       % start the  optimization routine
BVAR0                    = bvar_(y,lags,options);
pause;


%% Case 2: Estimation with mixed cherry picked-estimated Minnesota Prior
lags = 6; 
options.max_minn_hyper  = 1;       % start the  optimization routine
options.minn_prior_tau  = 10;      % set tau 
options.index_est       = [3 4];   % hyper-parameters to maximize
options.lb              = [0 0];   % lower bounds 
options.ub              = [20 20]; % upper bounds
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
% run BVAR with  optimal  hyperparameter selection
BVAR2                   = bvar_(y,lags,options);

% Plotting cholesky responses  to  Monetary  policy  shocks
% Define the IRF of Interest
indx_sho              = 4;   
% Change the order of the variables for the plot 
indx_var              = [4, 1, 2, 3];

% % IRFs to PLOT: compare IRFs obtained with different max hyperparam
mltple_irfs_to_plot_all(:,:,1,:) = BVAR0.ir_draws(indx_var,:,indx_sho,:);
mltple_irfs_to_plot_all(:,:,2,:) = BVAR1.ir_draws(indx_var,:,indx_sho,:);
mltple_irfs_to_plot_all(:,:,3,:) = BVAR2.ir_draws(indx_var,:,indx_sho,:);

% Customize the plot
% variables names for the plots
options.varnames      = {'1 year rate','IP','HICP','CORE INF'};  
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

%% ESSENCE: OPTIMALLY CHOSEN  AND  CHERRY  PICKED  CHOSEN  
%% HYPER-PARAMETER VALUES MAKE  LITTLE  DIFFERENCE  FOR  SHAPE  OF  IRFS.

%% 3.4 Use a presample as prior
clear options;
% run a VAR on presample data
presample = 50; % 8 years of presample
lags      = 6;
bvar1     = bvar_(y(1:presample,:),lags);
  

% use the VAR estimates to set the priors for the LP
options.priors.name        = 'Conjugate';
% posterior mean of the VAR AR coeff and constant
options.priors.Phi.mean    = mean(bvar1.Phi_draws,3);
% average variance of the AR coeff and constant
options.priors.Phi.cov     = diag(mean(var(bvar1.Phi_draws,0,3),2));
% posterior mean of the Covariance of the VAR residuals 
options.priors.Sigma.scale = mean(bvar1.Sigma_draws,3);
options.priors.Sigma.df    = size(bvar1.Phi_draws,1)-2;
options.K                  = 1000;

bvar2 = bvar_(y(presample+1:end,:),lags,options);




%% Case 4:  The Pandemic Priors
%           (with pandemic_on==1, the pandemic prior is active)
% example_pandemic_minn.m — VAR inference with Minnesota Prior vs Pandemic
% Priors(implementation of the Pandemic Priors of the Cascaldi-Garcia,2025)
%
% What this case does:
%   Estimates the same VAR twice with the BVAR_ toolbox: once as a
%   classic Minnesota Prior (no pandemic dummies), and once as a
%   Pandemic Priors (6 time dummies covering March-August
%   2020). Both runs share the same Minnesota hyperparameters, so the
%   only difference between the two estimations is the presence of the
%   pandemic dummies. 
%   First, we compare the log marginal data densities across the two specifications.
%   Subsequently, we plot the impulse responses to the EBP shock for each model
%   using plot_irfs_ to provide a direct visual comparison.
%   
%
% Reference:
%   Cascaldi-Garcia, D. (2025), "Pandemic Priors"
%   Paper: https://drive.google.com/file/d/1T0-q--zYZPRE_g1ijqL9NrKlh1X7Q16r/view
%   Website:  www.danilocascaldigarcia.com
%
%   The dataset used in this example (Data.xlsx) is taken directly from
%   the author's replication files
%% 4.1 Data preparation:
%   - Apply a log*100 transformation to the variables expressed in
%     levels (S&P 500, PCE, PCE Price Index, Employment, Industrial
%     Production), so that their coefficients/IRFs can be read as
%     approximate percentage changes; EBP, the Shadow Rate and the
%     Unemployment Rate are left in their original units (percentage
%     points), as in the original paper.
%   - Build a monthly date vector matching the sample, used later to
%     locate the start of the pandemic dummy window (March 2020).
data = readmatrix("Data.xlsx");
data = data(:,2:end);
Yraw = data;
log_vector = [0 1 0 1 1 1 1 0];
Yname = {'EBP','S&P 500','Shadow Rate','PCE','PCE Price Index','Employment','Ind. Production','Unemp. Rate'};
for ee = 1:size(log_vector,2)
    if log_vector(ee)==1; Yraw(:,ee) = log(Yraw(:,ee))*100; end
end
time_vec = datetime(1975,1,1):calmonths(1):datetime(2022,12,1);
nAR = 12;
covid_ind_F = find(datetime(2020,3,1)==time_vec);

%% 4.2 Log Marginal likelihood comparison

% ============ Common options ============
opts = struct();
opts.minn_prior_tau    = 5;
opts.minn_prior_decay  = 1;
opts.minn_prior_lambda = 0;
opts.minn_prior_mu     = 0.5;
opts.minn_prior_omega  = 1;
opts.K   = 2000;
opts.hor = 36;

% ============ Pandemic Priors ============
opts_pandemic = opts;
opts_pandemic.pandemic.start = covid_ind_F;    % first period covered by the pandemic dummies (March 2020)
opts_pandemic.pandemic.h     = 6;              % number of pandemic dummy periods (March-August 2020)
opts_pandemic.pandemic.phi   = 0.05;           % tightness of the pandemic dummies (small phi = uninformative,
% absorbs the anomaly; large phi = shrinks to zero, converging to the classic Minnesota Prior)

BVAR_pandemic = bvar_(Yraw, nAR, opts_pandemic);

% ============ Minnesota Prior  ============
opts_minn = opts;
BVAR_minn = bvar_(Yraw, nAR, opts_minn);

% ============ Log Marginal Likelihood Comparison ============
disp('log marginal data density')
disp(['  Minnesota Prior  = ' num2str(BVAR_minn.logmlike,'%6.2f')])
disp(['  Pandemic Priors  = ' num2str(BVAR_pandemic.logmlike,'%6.2f')])
pause;


%% 4.3 Comparison of EBP Shock Impact


% ============ EBP shock impact — Minnesota Prior ============
indx_sho = 1;              % EBP is the first variable -> the real EBP shock
indx_var = 1:size(Yraw,2); % all 8 variables

irfs_minn = BVAR_minn.ir_draws(indx_var,:,indx_sho,:);

clear plot_opts
plot_opts.varnames    = Yname;
plot_opts.shocksnames = {'EBP shock - Classic Minnesota'};
plot_irfs_(irfs_minn, plot_opts)
sgtitle('EBP shock - Classic Minnesota')

% ============ EBP shock impact — Pandemic Priors ============
irfs_pandemic = BVAR_pandemic.ir_draws(indx_var,:,indx_sho,:);

clear plot_opts
plot_opts.varnames    = Yname;
plot_opts.shocksnames = {'EBP shock - Pandemic Priors'};
plot_irfs_(irfs_pandemic, plot_opts)
sgtitle('EBP shock - Pandemic Priors')
pause;



%% Case 5:  Estimation with the Prior for the Long Run
%           (shrinks the combinations H*y, one tightness phi per row of H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Small phi imposes  a unit root on that combination, large phi is
% uninformative. Here the rows of H are the common trend Y+C+I and the two
% great ratios C-Y and I-Y: the variables are left free to trend, but not
% free to drift apart. The prior replaces the sum-of-coefficients and
% co-persistence dummies and nests them, since H = eye(ny) with
% phi = 1/minn_prior_mu gives back the minn_prior_mu block.
% Prior based on the paper "Priors for the Long Run" (2019) by D. Giannone,
% M. Lenza and G. Primiceri, JASA 114:526, 565-580
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

addpath ../../cmintools/
addpath ../../bvartools/

% load quarterly data from Giannone, Lenza and Primiceri (2019) data set
% US log real per capita GDP, consumption and investment, 1955Q1-2013Q1
load DataGLP
y    = [Y C I];   % the variables for the VAR
lags = 5;

%% 5.1 Estimation with a generic Minnesota prior
clear options
options.hor             = 60;   % long horizon: the prior acts on the long run
options.minn_prior_tau  = 3;    % overall tightness; decay, lambda, mu and
                                % omega stay at their default values
BVAR3                   = bvar_(y,lags,options);

%% 5.2 Estimation with the Minnesota and the long run prior
clear options
options.hor             = 60;
options.minn_prior_tau  = 3;    % same short run prior as in 5.1
options.priors.name     = 'PLR';
% the linear combinations the prior is elicited on
options.priors.PLR.H    = [ 1  1  1 ;     % Y+C+I, the common trend
                           -1  1  0 ;     % C-Y,   the consumption great ratio
                           -1  0  1 ];    % I-Y,   the investment great ratio
% one tightness per row of H; these are GLP's own estimates for this VAR
options.priors.PLR.phi  = [0.86; 0.56; 1.86];
BVAR4                   = bvar_(y,lags,options);

% Plotting cholesky responses to the first (GDP) shock
% Define the IRF of Interest
indx_sho              = 1;
% Order of the variables for the plot
indx_var              = [1, 2, 3];

% IRFs to PLOT: compare IRFs obtained with the two priors
mltple_irfs_to_plot_all(:,:,1,:) = BVAR3.ir_draws(indx_var,:,indx_sho,:);
mltple_irfs_to_plot_all(:,:,2,:) = BVAR4.ir_draws(indx_var,:,indx_sho,:);

% Customize the plot
clear options
% variables names for the plots
options.varnames      = {'GDP','Consumption','Investment'};
% name of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% name of the figure to save
options.saveas_strng  = 'minn_vs_longrun';
% name of the shock
options.shocksnames   = {'Minnesota','Minnesota + LR'};
% additional 90% HPD set
options.conf_sig_2    = 0.9;
% the plotting command
plot_all_irfs_(mltple_irfs_to_plot_all,options)

% The prior is elicited on H*y, so look at the responses of H*y itself: the
% two great ratios, with and without the long run prior
ratios_to_plot(1,:,1,:) = BVAR3.ir_draws(2,:,indx_sho,:) - BVAR3.ir_draws(1,:,indx_sho,:);
ratios_to_plot(2,:,1,:) = BVAR3.ir_draws(3,:,indx_sho,:) - BVAR3.ir_draws(1,:,indx_sho,:);
ratios_to_plot(1,:,2,:) = BVAR4.ir_draws(2,:,indx_sho,:) - BVAR4.ir_draws(1,:,indx_sho,:);
ratios_to_plot(2,:,2,:) = BVAR4.ir_draws(3,:,indx_sho,:) - BVAR4.ir_draws(1,:,indx_sho,:);

% Customize the plot
clear options
% variables names for the plots
options.varnames      = {'C - Y','I - Y'};
% name of the directory where the figure is saved
options.saveas_dir    = './irfs_plt';
% name of the figure to save
options.saveas_strng  = 'minn_vs_longrun_ratios';
% name of the shock
options.shocksnames   = {'Minnesota','Minnesota + LR'};
% additional 90% HPD set
options.conf_sig_2    = 0.9;
% the plotting command
plot_all_irfs_(ratios_to_plot,options)

%% ESSENCE: THE  LONG  RUN  PRIOR  BARELY  MOVES  THE  IMPACT  RESPONSES
%% BUT  DISCIPLINES  THE  LOW  FREQUENCIES:  THE  GREAT  RATIOS,  LEFT
%% FREE  TO  DRIFT  BY  THE  MINNESOTA  PRIOR,  ARE  PULLED  BACK
%% TOGETHER  AT  LONG  HORIZONS.
