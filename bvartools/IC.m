function [AIC, HQIC, BIC] = IC(S, E, T, K)
%IC Information criteria of a Gaussian VAR.
%
%   [AIC,HQIC,BIC] = IC(S,E,T,K)
%
%   Inputs
%     S   N x N   residual covariance
%     E   T x N   residuals
%     T   scalar  number of observations
%     K   scalar  number of estimated parameters
%
%   Outputs
%     AIC, HQIC, BIC   Akaike, Hannan-Quinn and Schwarz criteria, divided
%                      by T (lower is better)
%
%   See also OLS_REG, BVAR_.

%   Filippo Ferroni, 6/1/2015.

% E   = var.e_ols;
% S   = var.Sigma_ols; 

iS  = pinv(S);
N   = size(S,1);
llf = - (T * N / 2) * (1 + log(2 * pi)) - T / 2 * log(det(S));
llf = llf - 1 /2 * trace( iS * E' * E);

AIC = - 2 * llf / T + 2 * K / T;
% SIC = - 2 * llf / T + K * log(T) / T;
HQIC = - 2 * llf / T + 2 * K * log(log(T)) / T;
BIC = - 2 * llf / T + K * log(T) / T;
