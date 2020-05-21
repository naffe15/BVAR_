function [AIC, HQIC, BIC] = IC(var, T, K)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'IC' computes various information criteria

% Inputs:
% - llf, log marginal likelihood
% - T, sample size
% - K, number of regressors

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S   = var.Sigma_ols; iS  = inv(S);
E   = var.e_ols;
N   = size(S,1);
llf = - (T * N / 2) * (1 + log(2 * pi)) - T / 2 * log(det(S));
llf = llf - 1 /2 * trace( iS * E' * E);

AIC = - 2 * llf / T + 2 * K / T;
% SIC = - 2 * llf / T + K * log(T) / T;
HQIC = - 2 * llf / T + 2 * K * log(log(T)) / T;
BIC = - 2 * llf / T + K * log(T) / T;