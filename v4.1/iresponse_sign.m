function [ir,Omeg] = iresponse_sign(Phi,Sigma,hor,signrestriction,cont)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'iresponse_sign' computes the impulse response functions using sign
% restrictions on the endogenous variables
% Reference: Rubio-Ramirez, J. F., Waggoner, D. F. and Zha, T.: 2010,
% Structural Vector Autoregresions: Theory of Identification and Algorithms
% for Inference, Review of Economic Studies 77(2), 665–696.

% Inputs:
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - hor, horizon of the IRF
% - unit, 1 shock STD or 1 percent increase
% - signrestriction, cell array containing the restrictions. Sign
% restriction can be activated in the toolbox by setting the options
% options.signs{1} = 'y(a,b,c)>0'; 
% where a, b and c are integer. The syntax means that shock c has a
% positive impact on the a-th variable at horizon b.  

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
tol     = 0;
favar   = 1;

if nargin < 5
    cont  = eye(n);
    favar = 0;
end


while d==0 && tol < 30000
    %     e = randn(n,1);
    %     Omega = (e./norm(e));
    %
    G     = randn(m);
    [Q,R] = qr(G);
    % normalize to positive entry in the diagonal
    In    = diag(sign(diag(R)));
    Omega = Q  * In;
    %    u   = Ao * Q * In;
    y = iresponse(Phi,Sigma,hor,Omega);
    if favar == 1 
        for jj  = 1 : size(y,3) 
            yy(:,:,jj) = cont * y(:,:,jj);
        end
        clear y; y = yy;
    end
    count   = 0;
    for ii = 1 : size(signrestriction,2),
        tmp = eval(signrestriction{ii});
        count = count + tmp;
    end
    if count == size(signrestriction,2), % if all signs are verified stop
        d=1;
        Omeg  = Omega;
        ir    = y;
        %      if tmp==1, % if all signs are verified stop
        %         d=1;
        %         Omeg  = Omega;
        %         ir    = y;
    else
        d=0; tol = tol + 1;
    end
end

if d==0
%     ir = iresponse(Phi,Sigma,hor,Omega);
% else
    warning('I could not find a rotation satisfying the restrictions.')
end
end
