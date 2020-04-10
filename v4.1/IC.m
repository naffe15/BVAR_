function [AIC, SIC,HQIC, BIC] = IC(llf, T, K)
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


if (nargin == 1) && isstruct(llf)
    T   = llf.T;
    K   = llf.K;
    llf = llf.llf;
end

AIC = - 2 * llf / T + 2 * K / T;
SIC = - 2 * llf / T + K * log(T) / T;
HQIC = - 2 * llf / T + 2 * K * log(log(T)) / T;
BIC = - 2 * llf / T + K * log(T);