function [ir,Q]=iresponse_longrun(alpha,Sigma,hor,lags)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse_longrun' computes the impulse response functions to a long run shock
% ordered first 

% Inputs:
% - alpha, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - hor, horizon of the IRF
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

N            = size(Sigma,1);
ir          = zeros(N,hor,N); % variables, horizon, shock

% companion
F       = [alpha(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
Inp     = eye(size(F));
% long run sum of VMA 
C1             = Inp - F;
C1             = C1 \ Inp; % = inv(C1)*Inp;
C1             = C1(1 : N, 1 : N);
PSI1           = chol(C1 * Sigma * C1')';
% IMPACT MATRIX
Q              = (C1 \ eye(size(C1))) * PSI1;
% normalize to positive entry in the element Q(1,1)
Q(1,1)     = Q(1,1) * sign(Q(1,1));
Fk         = eye(N * lags);

for k=1:hor
    PHI         = Fk * G * Q;
    ir(1:N,k,:) = G'* PHI;
    Fk          = F * Fk;
end

    