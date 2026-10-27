function [ir,Omeg] = iresponse_sign_hmoments(errors,Phi,Sigma,hor,signrestriction,hmoments,f,cont)
%IRESPONSE_SIGN_HMOMENTS Responses identified with signs and higher-moment restrictions.
%
%   [ir,Omeg] = IRESPONSE_SIGN_HMOMENTS(errors,Phi,Sigma,hor,signrestriction,hmoments,f)
%   [ir,Omeg] = IRESPONSE_SIGN_HMOMENTS(...,cont)
%
%   Draws random rotations (up to 30000) until the responses satisfy the
%   sign restrictions (as drawn or with signs flipped) and the implied
%   structural shocks satisfy the restrictions on their skewness and
%   kurtosis.
%
%   Inputs
%     errors           T x ny     reduced-form residuals
%     Phi              (ny*lags+nx) x ny   VAR coefficients (see IRESPONSE)
%     Sigma            ny x ny    covariance of the reduced-form residuals
%     hor              scalar     number of horizons
%     signrestriction  cell array of strings 'y(a,b,c)>0' (options.signs)
%     hmoments         cell array of strings 'hm(a,b,c) > threshold'
%                      (options.hmoments; see CHECKMRESTRICTIONS)
%     f                ny x 4 x 2 from HMOMENTS2MATRIX
%     cont             nobs x ny  FAVAR only: loadings of the observables
%
%   Outputs
%     ir     ny x hor x ny   responses (nobs x hor x ny with cont), NaN if
%                            no rotation is found
%     Omeg   ny x ny         accepted rotation of chol(Sigma,'lower')
%
%   See also CHECKMRESTRICTIONS, HMOMENTS2MATRIX, IRESPONSE_SIGN.
%
%   Reference: Andrade, P., F. Ferroni and L. Melosi (2024); see Section
%   3.2.9 of the Hitchhiker's Guide.

%   Filippo Ferroni.

[m,n]   = size(Sigma);
d       = 0;
% d0      = 0;
tol     = 0;
favar   = 1;
if nargin < 8
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

A       = chol(Sigma,'lower');
v       = zeros(size(errors));
hor0    = maxhorsign(signrestriction) ;


while d==0 && tol < 30000
    % generate a random orthonormal matrix
    Omega = generateQ(m);
    % compute IRF
    y = iresponse(Phi,Sigma,hor0,Omega);
    % uncompress the factors (if favar)
    if favar == 1
        if do_pagemtimes == 0
            for jj  = 1 : size(y,3)
                yy(:,1:hor0,jj) = cont * y(:,1:hor0,jj);
            end
        else
            lam  = repmat(cont,1,1,n);
            yy   = pagemtimes(lam,y);
        end
        clear y; y = yy;
    end
    % check sign restrictions
    fsign = 1;
    %d0 = checkrestrictions(signrestriction,y);
    [d0, fsign] = checkrestrictions2(signrestriction,y);
    if d0 ==1
        v = errors / ( (fsign * Omega)' * A');  % structural innovations
        % check higher-moments restrictions
        d = checkmrestrictions(hmoments,v, f);
    end
    if d==1 % stop 
%         Omeg  = Omega;
        Omeg  = fsign * Omega;
        y     = iresponse(Phi,Sigma,hor,Omeg);
        ir    = y;
        if favar == 1   % map the responses into the observables
            for jj  = 1 : size(y,3)
                yy(:,1:hor,jj) = cont * ir(:,1:hor,jj);
            end
            clear ir; ir = yy;
        end
    else
        tol = tol + 1;
    end
end

if d==0
    warning('I could not find a rotation satisfying the restrictions.')
end
end
