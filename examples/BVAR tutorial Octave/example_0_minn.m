%% BVAR tutorial: Inference with Minnesota Prior
% Author:   Filippo Ferroni
% Date:     27/02/2020

clear all
close all
clc
addpath ../../cmintools/
addpath ../../v4.2/
pkg load statistics

%% %=========================================================================
%%% INFERENCE %%%
%%=========================================================================

% load the data
load('../BVAR tutorial/Data.mat')
y= [IPI HICP CORE Euribor1Y M3 EXRATE];% collect the variables used in the VAR

%% Ex1/ Minnesota Prior

lags = 6; 
options.max_minn_hyper  = 1;
options.minn_prior_tau  = 10;      % set tau 
options.index_est       = [3 4];   % define the hyper-parameters over which to maximize
options.lb              = [0 0];   % sets the lower bounds 
options.ub              = [20 20]; % sets the upper bounds
options.max_compute     = 3;       % optimization  by Matlab Simplex
BVAR                    = bvar(y,lags,options);

%% Ex2/ Minnesota Prior

clear options
lags = 6; 
% setting the default values for the hyperparameters
hyperpara(1)    = 3;		  % tau
hyperpara(2)    = 0.5;		  % decay
hyperpara(3)    = 5;		  % lambda
hyperpara(4)    = 2;		  % mu
hyperpara(5)    = 2;		  % omega
% setting the options
options.index_est	   = 1:1;      % hyper-parameter over which maximize
options.max_compute    = 3;      % maximize  using Matlab fmincon function
options.lb             = [0.05]; % Lower bound
[postmode,logmlike,HH] = bvar_max_hyper(hyperpara,y,lags,options); 


%% Ex3

hyperpara(1)            = postmode(1); % use as starting value previous mode
options.index_est       = 1:3; % set hyper-parameters over which maximize
options.lb              = [0.05 0.05 0.05]; % Lower bounds
options.ub              = [50 50 50];       % Upper bounds
[postmode,log_dnsty,HH] = bvar_max_hyper(hyperpara,y,lags,options);
