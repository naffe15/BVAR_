function XLag = lagX(X,lags)
%lagX Create matrix of lagged time series, 
% 1st dimension time, 
% 2nd dimention variables
%
%if size(X,2) > size(X,1)
%   X = X'; % Ensure a column vector
%end

missingValue = NaN;  % Assign default missing value
Lags = length(lags); % Number of lags to apply to each time series
[T,ny] = size(X);
XLag = missingValue(ones(T,ny*Lags)); % Preallocate

for c = 1:Lags

    L       = lags(c);
    columns = (ny*(c-1)+1):c*ny; % Columns to fill, this lag

    if L > 0 % Time delays
       XLag((L + 1):end,columns) = X(1:(end - L), :);

    elseif L < 0 % Time leads
       XLag(1:(end + L),columns) = X((1 - L):end, :);

    else % No shifts
       XLag(:,columns) = X;

    end

end
