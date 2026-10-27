function [C] = connectedness(Phi,Sigma,p,Omega)
%CONNECTEDNESS Spillover (connectedness) measures of Diebold and Yilmaz.
%
%   C = CONNECTEDNESS(Phi,Sigma,p)
%   C = CONNECTEDNESS(Phi,Sigma,p,Omega)
%
%   Based on the p-step forecast error variance decomposition: generalised
%   (Pesaran-Shin) when Omega is omitted, otherwise the one of the
%   structural shocks with impact matrix Omega.
%
%   Inputs
%     Phi     (ny*lags) x ny   autoregressive coefficients only, lag-1 block
%                              first (no constant)
%     Sigma   ny x ny          residual covariance
%     p       scalar           horizon (options.nethor in BVAR_)
%     Omega   ny x ny          impact matrix, e.g. chol(Sigma,'lower')*Q
%
%   Output
%     C.Index           scalar   total connectedness, percent
%     C.FromAllToUnit   ny x 1   received by each variable from the others
%     C.FromUnitToAll   ny x 1   transmitted by each variable to the others
%     C.Net             ny x 1   transmitted minus received
%     C.theta           ny x ny  variance shares before normalisation
%
%   See also VAR2MA, BVAR_.
%
%   Reference: Diebold, F. X. and K. Yilmaz (2012), "Better to Give than to
%   Receive: Predictive Directional Measurement of Volatility Spillovers",
%   International Journal of Forecasting 28(1), 57-66.

%   Filippo Ferroni, 5/1/2021.

Pesaran_Shin = 0;
% dimensions
[m, ny] = size(Phi);    
if rem(m, ny)==0
    lags = m/ny;
else
    lags = floor((m-1)/ny);
end
% if no restriction are imposed use generalized
if nargin<4
    % identification matrix, if nothing is declared use the Pesaran-Shin
    % generalized IRF
    Pesaran_Shin  = 1; 
    Omega = Sigma;        
end

% preaallocation
theta = zeros(ny);
% remove the constant/time trend from the Phi matrix
Phi0 = Phi(1:lags*ny,:);
% compute the MA representation
Psi = var2ma(Phi0,p);
% Equation (1) in Diebold and Yilmaz (2012)
% Construct the Denominator
DD = 0;
for hh = 0 : p-1
    DD = DD + Psi(:,:,hh+1) * Sigma * Psi(:,:,hh+1)';
    % DD = DD + Psi(:,:,hh+1) * Omega * Psi(:,:,hh+1)';
end
% Construct the Numerator and theta
for ss = 1 : ny
    for j= 1 : ny
        NN = 0;
        for hh = 0 : p-1
            %AA = Psi(:,:,hh+1)*Sigma;
            AA = Psi(:,:,hh+1) * Omega;
            NN = NN + AA(ss,j)^2;
        end
        theta(ss,j) = NN / DD(ss,ss);
        if Pesaran_Shin  == 1 
            theta(ss,j) = Sigma(j,j)^(-1/2) * theta(ss,j);            
        end
    end
end

% normalize by the sum of the rows
Theta = theta ./ repmat(sum(theta,2),1,ny);

if max(max(abs(sum(Theta,2) - ones(ny,1))))>1e-8
    error('something went wrong')
end

% remove the main diagonal
% Theta0          = triu(Theta)-tril(Theta);
Theta0          = Theta - diag(diag((Theta)));
C.Index         = sum(sum(Theta0)) / ny * 100; % overall index of connectedness
C.FromAllToUnit = sum(Theta0,2) / (ny-1) * 100;    % connectedness from all to unit i 
C.FromUnitToAll = sum(Theta0,1)' / (ny-1) * 100;   % connectedness from unit i to all 
C.Net           = C.FromUnitToAll - C.FromAllToUnit;   % connectedness from unit i to all 
C.theta         = theta; 
