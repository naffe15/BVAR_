function [ir,Omeg] = iresponse_sign_narrative(errors,Phi,Sigma,hor,signrestriction,narrative,cont)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse_sign_narrative' computes the impulse response functions using sign
% restrictions on the endogenous variables
% Reference:

% Inputs:
% - e, VAR reduced form errors (T x n)
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - hor, horizon of the IRF
% - unit, 1 shock STD or 1 percent increase
% - signrestriction, cell array containing the restrictions. Sign
% restriction can be activated in the toolbox by setting the options
% options.signs{1} = 'y(a,b,c)>0';
% where a, b and c are integer. The syntax means that shock c has a
% positive impact on the a-th variable at horizon b.
% - narrative, cell array containing the restrictions. Sign
% restriction can be activated in the toolbox by setting the options
% options.narrative{1} = 'v(a,b)>0';
% where m can be a vector of scalar and n is a scalar. The syntax means
% that shock n has to be positive .

% Output:
% - ir contains the IRF
% 1st dimension:   variable
% 2st dimension:   horizon
% 3st dimension:   shock

% Filippo Ferroni, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
                yy(:,:,jj) = cont * y(:,:,jj);
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
