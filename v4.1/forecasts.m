function [frcst_no_shock,frcst_with_shocks]=forecasts(forecast_data,Phi,Sigma,fhor,lags,EPS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'forecasts' computes out-of-sample forecasts

% Inputs:
% - forecast_data, last data
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - fhor, horizon of the out-of-sample forecasts
% - lags
% - EPS, a particular shock realization

% Output:
% - frcst_no_shock
% 1st dimension:   horizon 
% 2nd dimension:   variable 

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ny = size(Sigma,1);
if nargin < 6
    shock_given = 0;
else
    shock_given = 1;
end

Sigma_lower_chol = chol(Sigma)';

% preallocating memory
frcst_no_shock     = nan(fhor,ny);
frcst_with_shocks  = nan(fhor,ny);

lags_data = forecast_data.initval;
for t = 1 : fhor
    X = [ reshape(flip(lags_data, 1)', 1, ny*lags) forecast_data.xdata(t, :) ];
    y = X * Phi;
    lags_data(1:end-1,:) = lags_data(2:end, :);
    lags_data(end,:) = y;
    frcst_no_shock(t, :) = y;
end

% With shocks
lags_data = forecast_data.initval;
for t = 1 : fhor
    X = [ reshape(flip(lags_data, 1)', 1, ny*lags) forecast_data.xdata(t, :) ];
    shock = (Sigma_lower_chol * randn(ny, 1))';
    if shock_given == 1
        shock  = EPS(t,:);
    end
    y = X * Phi + shock;
    lags_data(1:end-1,:) = lags_data(2:end, :);
    lags_data(end,:) = y;
    frcst_with_shocks(t, :) = y;
end

