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

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[m,n]   = size(Sigma);
ir      = nan(n,hor,n);
Omeg    = nan(n);
d       = 0;
d0      = 0;
tol     = 0;
favar   = 1;
if nargin < 7
    cont  = eye(n);
    favar = 0;
end


A       = chol(Sigma,'lower');
v       = zeros(size(errors));

while d==0 && tol < 30000
    % generate a random orthonormal matrix
    Omega = generateQ(m);
    % compute IRF
    y = iresponse(Phi,Sigma,hor,Omega);
    % uncompress the factors (if favar)
    if favar == 1
        for jj  = 1 : size(y,3)
            yy(:,:,jj) = cont * y(:,:,jj);
        end
        clear y; y = yy;
    end
    % check sign restrictions
    d0 = checkrestrictions(signrestriction,y);
    if d0 ==1
        v = errors / ( Omega' * A');  % structural innovations
        % check narrative restrictions
        d = checkrestrictions(narrative,y,v);
    end
    if d==1 % stop 
        Omeg  = Omega;
        ir    = y;
    else
        tol = tol + 1;
    end
end

if d==0
    warning('I could not find a rotation satisfying the restrictions.')
end
end
