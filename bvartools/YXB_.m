function [YYact,XXact] = YXB_(YY,lags,constant_timetrend)
%YXB_ Builds the matrices of the regression Y = X*B + E for a VAR.
%
%   [YYact,XXact] = YXB_(YY,lags)
%   [YYact,XXact] = YXB_(YY,lags,constant_timetrend)
%
%   Inputs
%     YY                  T x nv   data
%     lags                scalar   number of lags
%     constant_timetrend  1 x 2    [constant timetrend], 1 to include each
%                                  (default [1 0])
%
%   Outputs
%     YYact   (T-lags) x nv              left-hand side
%     XXact   (T-lags) x (nv*lags+d)     lags (lag-1 block first), then the
%                                        constant and the trend
%
%   No dummy observations are added.
%
%   See also RFVAR3, FAST_OLS_.

%   Filippo Ferroni, 6/1/2015.

if nargin<3
    constant = 1;
    timetrend = 0;
else 
    constant  = constant_timetrend(1);
    timetrend = constant_timetrend(2);
end

nlags_   = lags;                 % number of lags   */
T0       = lags;                 % size of pre-sample */

nv      = size(YY,2);            %* number of variables */
nobs    = size(YY,1)-T0;         %* number of observations */

% Actual observations

YYact = YY(T0+1:T0+nobs,:);
XXact = zeros(nobs,nv*nlags_);
i = 1;

while (i <= nlags_)
    XXact(:,(i-1)*nv+1:i*nv) = YY(T0-(i-1):T0+nobs-i,:);
    i = i+1;
end

if constant
    % last column of XXact = constant
    XXact = [XXact ones(nobs,1)];
end

if timetrend
    % last column of XXact = constant
    XXact = [XXact (1:nobs)'];
end

    
