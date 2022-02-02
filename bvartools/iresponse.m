function [ir]=iresponse(alpha,Sigma,hor,Omega,unit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse' computes the impulse response functions

% Inputs:
% - alpha, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - hor, horizon of the IRF
% - Omega, a particular orthonormal rotation
% - unit, 1 shock STD or 1 percent increase

% Output:
% - ir contains the IRF 
% 1st dimension:   variable 
% 2st dimension:   horizon 
% 3st dimension:   shock

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N           = size(Sigma,1);
% ir          = zeros(N,hor,N); % variables, horizon, shock
ir          = zeros(N,hor,size(Omega,2)); % variables, horizon, shock/exogenous variable
[m , n]     = size(alpha);
if rem(m, n)==0
    lags = m/n;
else
    lags = floor((m-1)/n);
end
try 
    [Q]         = chol(Sigma,'lower');
catch
    iS               = pinv(Sigma);
    [L, D, P]        = ldl(iS);
    S_inv_upper_chol = sqrt(D) * L' * P';
    Q = S_inv_upper_chol';
end


% units 
if nargin < 5
    unit = eye(N);
else
    unit = inv(diag(diag(Q)));
end
    
% 1 standard deviation increase (if unit = eye(N))
% 100 basis point increase (else)
Q = Q*unit;

% companion form 
F       = [alpha(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
Fk      = eye(N * lags);
% compute IRFs
for k=1:hor
    PHI         = Fk * G * Q * Omega;
    ir(1:N,k,:) = G' * PHI;
    Fk          = F * Fk;
end
    