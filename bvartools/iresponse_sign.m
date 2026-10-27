function [ir,Omeg] = iresponse_sign(Phi,Sigma,hor,signrestriction,cont)
%IRESPONSE_SIGN Impulse responses identified with sign restrictions.
%
%   [ir,Omeg] = IRESPONSE_SIGN(Phi,Sigma,hor,signrestriction)
%   [ir,Omeg] = IRESPONSE_SIGN(Phi,Sigma,hor,signrestriction,cont)
%
%   Draws random rotations (up to 30000) until one satisfies all the
%   restrictions, either as drawn or with all signs flipped.
%
%   Inputs
%     Phi              (ny*lags+nx) x ny   VAR coefficients (see IRESPONSE)
%     Sigma            ny x ny             covariance of the reduced-form
%                                          residuals
%     hor              scalar              number of horizons
%     signrestriction  cell array of strings, options.signs in BVAR_, e.g.
%                      'y(a,b,c)>0': shock c raises variable a at horizon
%                      b (b = 1 is the impact)
%     cont             nobs x ny           FAVAR only: loadings mapping the
%                                          VAR variables into the observables
%
%   Outputs
%     ir     ny x hor x ny   responses (nobs x hor x ny with cont), NaN if
%                            no rotation is found
%     Omeg   ny x ny         accepted rotation of chol(Sigma,'lower')
%
%   See also IRESPONSE_SIGN_NARRATIVE, IRESPONSE_SET, CHECKRESTRICTIONS2,
%   GENERATEQ.
%
%   Reference: Rubio-Ramirez, J. F., D. F. Waggoner and T. Zha (2010),
%   "Structural Vector Autoregressions: Theory of Identification and
%   Algorithms for Inference", Review of Economic Studies 77(2), 665-696.

%   Filippo Ferroni.

[m,n]   = size(Sigma);
d       = 0;
tol     = 0;
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
hor0 = maxhorsign(signrestriction) ;

    
while d==0 && tol < 30000
    % generate a random orthonormal matrix
    Omega = generateQ(m);
    % compute IRF
%     y = iresponse(Phi,Sigma,hor,Omega);
    y = iresponse(Phi,Sigma,hor0,Omega);
    % uncompress the factors (if favar)
    if favar == 1   
        if do_pagemtimes ==0
            for jj  = 1 : size(y,3)
                yy(:,1:hor0,jj) = cont * y(:,1:hor0,jj);
            end
        else
            lam  = repmat(cont,1,1,n);
            yy   = pagemtimes(lam,y);
        end
        clear y; y = yy;        
    end
    % check restrictions
    % fsign = 1;
    % d = checkrestrictions(signrestriction,y);
    [d,fsign] = checkrestrictions2(signrestriction,y);
    if d==1
        Omeg  = fsign * Omega;
        y     = iresponse(Phi,Sigma,hor,Omega);
        ir    = fsign * y;
        if favar == 1
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
