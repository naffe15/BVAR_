function MA = var2ma(AR,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% var2ma coverts the VAR into a MA
% AR only contains the AR coefficients and does not contain constant, time
% trends or exogenous variables. 
% Filippo Ferroni, 5/1/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ny   = size(AR,2);
lags = size(AR,1)/ny;
MA        = nan(ny,ny,p);
MA(:,:,1) = eye(size(AR,2));
MA(:,:,2) = AR(1:ny, 1:ny);

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
        tmp    = tmp + MA(:, :, spanma) * AR(spanar, :) ;
    end
    MA(:,:,jj) = tmp;
end

