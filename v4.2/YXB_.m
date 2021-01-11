function [YYact,XXact] = YXB_(YY,lags,constant_timetrend)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'YXB_' organizes the data in the form
% of Y = XB+E
% NO dummy observations
% constant and trends are at the end

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

    
