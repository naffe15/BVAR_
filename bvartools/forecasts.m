function [frcst_no_shock,frcst_with_shocks]=forecasts(forecast_data,Phi,Sigma,fhor,lags,EPS,skip_)
%FORECASTS Out-of-sample forecasts of a VAR, with and without shocks.
%
%   [frcst_no_shock,frcst_with_shocks] = FORECASTS(forecast_data,Phi,Sigma,fhor,lags)
%   [...] = FORECASTS(forecast_data,Phi,Sigma,fhor,lags,EPS)
%   [...] = FORECASTS(forecast_data,Phi,Sigma,fhor,lags,EPS,skip_)
%
%   Inputs
%     forecast_data  struct with fields
%       .initval     lags x ny x nunits    last lags observations, oldest first
%       .xdata       fhor x nx             exogenous regressors over the forecast
%                                          horizon (constant included)
%     Phi            (ny*lags+nx) x ny     VAR coefficients, lag-1 block first,
%                                          exogenous rows last
%     Sigma          ny x ny               covariance of the reduced-form residuals
%     fhor           scalar                forecast horizon
%     lags           scalar                number of lags
%     EPS            fhor x ny             residual path for the with-shocks
%                                          forecast (same for every unit). If
%                                          omitted, residuals are drawn from
%                                          N(0,Sigma).
%     skip_          0 (default) or 1      1 skips the no-shock forecast, which
%                                          is then returned as NaN
%
%   Outputs
%     frcst_no_shock     fhor x ny x nunits   path with future residuals set to 0
%     frcst_with_shocks  fhor x ny x nunits   path with residuals EPS or one
%                                             random draw
%
%   See also CFORECASTS, CFORECASTS2, PLOT_FRCST_, BVAR_.

%   Filippo Ferroni, 6/1/2015.

ny = size(Sigma,1);
nunits = size(forecast_data.initval,3);

if nargin < 6
    shock_given = 0;
else
    shock_given = 1;
end
if nargin < 7
    skip_ = 0;
end

Sigma_lower_chol = chol(Sigma)';

% preallocating memory
frcst_no_shock     = nan(fhor,ny,nunits);
frcst_with_shocks  = nan(fhor,ny,nunits);

for uu  = 1 : nunits   
    % no shocks
    lags_data = forecast_data.initval(:,:,uu);
    if skip_ == 0
        for t = 1 : fhor
            X = [ reshape(flip(lags_data, 1)', 1, ny*lags) forecast_data.xdata(t, :) ];
            y = X * Phi;
            lags_data(1:end-1,:) = lags_data(2:end, :);
            lags_data(end,:) = y;
            frcst_no_shock(t, :, uu) = y;
        end
    end

    % With shocks
    lags_data = forecast_data.initval(:,:,uu);
    for t = 1 : fhor
        X = [ reshape(flip(lags_data, 1)', 1, ny*lags) forecast_data.xdata(t, :) ];
        if shock_given == 1
            shock  = EPS(t,:);
        else
            shock = (Sigma_lower_chol * randn(ny, 1))';
        end
        y = X * Phi + shock;
        lags_data(1:end-1,:) = lags_data(2:end, :);
        lags_data(end,:) = y;
        frcst_with_shocks(t, :, uu) = y;
    end

end