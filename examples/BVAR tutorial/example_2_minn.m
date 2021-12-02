% example_2_min_new.m Inference with Minnesota Prior
% Author:   Filippo Ferroni and  Fabio Canova
% Date:     27/02/2020, revised  14/10/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute  posterior estimate with calibrated Minnesota prior and  with
% optimally chosen  Minnesota prior.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
warning off;
close all; clc;

addpath ../../cmintools/
addpath ../../bvartools/


% load the Euro data
load Data
y= [IPI HICP CORE Euribor1Y M3 EXRATE];% the variables for the VAR

%% Case 1: Estimation with calibrated/ estimated Minnesota Prior

lags = 6; 
options.max_minn_hyper  = 1;       % start the  optimization routine
options.minn_prior_tau  = 10;      % set tau 
options.index_est       = [3 4];   % hyper-parameters to maximize
options.lb              = [0 0];   % the lower bounds 
options.ub              = [20 20]; % the upper bounds
options.max_compute     = 7;       % optimization  by Matlab Simplex
BVAR                    = bvar_(y,lags,options);
pause;

%% Case 2:  Optimization of one dimension of the Minnesota Prior
%           (without  posterior draws to speed up computations) 
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

%% Case 3: Take  previous optimal value  and compute optimal values for two
%             parameters (without  posterior  draws)

hyperpara(1)            = postmode(1); % use as starting value previous mode
options.index_est       = 1:3;         % set hyper-parameters over which maximize
options.lb              = [0.1 0.1 0.1]; % Lower bounds
options.ub              = [10 10 10];    % Upper bounds
[postmode,log_dnsty,~] = bvar_max_hyper(hyperpara,y,lags,options);

