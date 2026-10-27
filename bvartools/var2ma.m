function MA = var2ma(AR,p)
%VAR2MA Moving-average coefficients of a VAR.
%
%   MA = VAR2MA(AR,p)
%
%   Inputs
%     AR   (ny*lags) x ny   autoregressive coefficients only, lag-1 block
%                           first (no constant, trend or exogenous rows)
%     p    scalar           number of MA matrices (>= 2); p can exceed lags
%
%   Output
%     MA   ny x ny x p      MA(:,:,j) is the coefficient on the residual
%                           j-1 periods back; MA(:,:,1) = eye(ny)
%
%   See also CONNECTEDNESS, FE_, IRESPONSE.

%   Filippo Ferroni, 5/1/2021.

ny   = size(AR,2);
lags = size(AR,1)/ny;
MA        = nan(ny,ny,p);
MA(:,:,1) = eye(size(AR,2));
MA(:,:,2) = AR(1:ny, 1:ny)';

for jj = 3 : p
    tmp = 0;    
    if jj <= lags
        Kstp = jj;
    else
        Kstp = lags+1;
    end
    for ll = 2 : Kstp
        spanma = jj - ll + 1;
        spanar = (ll-2)*ny+1 : (ll-2)*ny + ny;
        tmp    = tmp + MA(:, :, spanma) * AR(spanar, :)' ;
    end
    MA(:,:,jj) = tmp;
end

