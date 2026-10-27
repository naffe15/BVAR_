function [ir,Omeg] = iresponse_sign_narrative(errors,Phi,Sigma,hor,signrestriction,narrative,cont)
%IRESPONSE_SIGN_NARRATIVE Responses identified with sign and narrative restrictions.
%
%   [ir,Omeg] = IRESPONSE_SIGN_NARRATIVE(errors,Phi,Sigma,hor,signrestriction,narrative)
%   [ir,Omeg] = IRESPONSE_SIGN_NARRATIVE(errors,Phi,Sigma,hor,signrestriction,narrative,cont)
%
%   Draws random rotations (up to 30000) until one satisfies the sign
%   restrictions on the responses and the narrative restrictions on the
%   implied structural shocks.
%
%   Inputs
%     errors           T x ny     reduced-form residuals
%     Phi              (ny*lags+nx) x ny   VAR coefficients (see IRESPONSE)
%     Sigma            ny x ny    covariance of the reduced-form residuals
%     hor              scalar     number of horizons
%     signrestriction  cell array of strings 'y(a,b,c)>0' (options.signs)
%     narrative        cell array of strings in v (options.narrative), e.g.
%                      'v(tau,n)>0': shock n is positive in period tau (a
%                      row of errors); 'sum(v(tau0:tau1,n))>0' is allowed
%     cont             nobs x ny  FAVAR only: loadings of the observables
%
%   Outputs
%     ir     ny x hor x ny   responses (nobs x hor x ny with cont), NaN if
%                            no rotation is found
%     Omeg   ny x ny         accepted rotation of chol(Sigma,'lower')
%
%   No importance weights are computed; IRESPONSE_SET returns them.
%
%   See also IRESPONSE_SIGN, IRESPONSE_SET, CHECKRESTRICTIONS.
%
%   Reference: Antolin-Diaz, J. and J. F. Rubio-Ramirez (2018), "Narrative
%   Sign Restrictions for SVARs", American Economic Review 108(10),
%   2802-2829.

%   Filippo Ferroni.

[m,n]   = size(Sigma);
d       = 0;
d0      = 0;
tol     = 0;
favar   = 1;
if nargin < 7
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
    %y = iresponse_(Phi,A,hor,Omega);
    %y = iresponse(Phi,Sigma,hor,Omega);
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
    d0 = checkrestrictions(signrestriction,y);
    %[d0, fsign] = checkrestrictions2(signrestriction,y);
    if d0 ==1
        v = errors / ( (fsign * Omega)' * A');  % structural innovations
        % check narrative restrictions
        d = checkrestrictions(narrative,[],v);
    end
    if d==1 % stop 
        Omeg  = Omega;
        y = iresponse(Phi,Sigma,hor,Omega);
        ir    = y;
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
