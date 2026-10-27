function [ir,Omeg] = iresponse_heterosked(Phi,errors,hor,delta,cont)
%IRESPONSE_HETEROSKED Responses identified through heteroskedasticity (two regimes).
%
%   [ir,Omeg] = IRESPONSE_HETEROSKED(Phi,errors,hor,delta)
%   [ir,Omeg] = IRESPONSE_HETEROSKED(Phi,errors,hor,delta,cont)
%
%   The rotation comes from the singular value decomposition of
%   inv(S0)*S1, where S0 and S1 are the residual covariances in the two
%   regimes. Used by BVAR_ with options.heterosked_regimes.
%
%   Inputs
%     Phi      (ny*lags+nx) x ny   VAR coefficients (see IRESPONSE)
%     errors   T x ny    reduced-form residuals
%     hor      scalar    number of horizons
%     delta    T x 1     regime indicator: 0 first regime, 1 second regime
%     cont     nobs x ny FAVAR only: loadings of the observables
%
%   Outputs
%     ir       ny x hor x ny   responses (nobs x hor x ny with cont)
%     Omeg     ny x ny         rotation
%
%   See also IRESPONSE.
%
%   Reference: Sims, C. A. (2020), "SVAR Identification Through
%   Heteroskedasticity with Misspecified Regimes", Princeton University.

%   Filippo Ferroni.

[~,n]   = size(Phi);
favar   = 1;
if nargin < 5
    cont  = eye(n);
    favar = 0;
end
Omeg    = nan(n);
if favar == 0
    ir      = nan(n,hor,n);
else
    ir      = nan(size(cont,1),hor,n);
    yy      = nan(size(cont,1),hor,n);   
end
do_pagemtimes = 0;
if exist('pagemtimes','builtin') == 5
    do_pagemtimes = 1;
end

% compute the impact matrix 
S0      = cov(errors((delta==0),:));
S1      = cov(errors((delta==1),:));
X       = S0\S1;  %iS0S1 = inv(S0')*S1;
[~,S,V] = svd(X);
Omega   = V;
% compute impulse responses
y = iresponse(Phi,S,hor,Omega);
% uncompress the factors (if favar)
if favar == 1
    if do_pagemtimes ==0
        for jj  = 1 : size(y,3)
            yy(:,:,jj) = cont * y(:,:,jj);
        end
    else
        lam  = repmat(cont,1,1,n);
        yy   = pagemtimes(lam,y);
    end
    clear y; y = yy;
end
ir = y;
Omeg = Omega;