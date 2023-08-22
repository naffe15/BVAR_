function [ir,Omeg] = iresponse_sign_hmoments_eig(errors,Phi,Sigma,hor,moment)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

% Inputs:
% - e, VAR reduced form errors (T x n)
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - hor, horizon of the IRF

% Output:
% - ir contains the IRF
% 1st dimension:   variable
% 2st dimension:   horizon
% 3st dimension:   shock

% Filippo Ferroni, 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n]   = size(Sigma);
if nargin < 5
    moment  = 4;    
end

A       = chol(Sigma,'lower');
%  varols.u * BVAR.iSigChol';
u     = errors * A';

if moment == 4
    % eigenvalue decompostion of fourth moments
    Kz    = eye(n^2) + commutationmatrix(n) + reshape(eye(n),n*n,1)*reshape(eye(n),n*n,1)';
    Kappa = fourthmom(u,1) - Kz;
    [P4, ~, ~] = svd(Kappa);
    OO         = nan(n);
    for ni  = 1 : n
        OO(:,ni) = sign(P4(1,ni)) * P4(1:n,ni) ./ sqrt(abs(P4(1,ni)));
    end
elseif moment ==3
    % eigenvalue decompostion of third moments
    S = thirdmom(u,1);
    [OO, ~, ~] = svd(S*S');    
end

Omeg = OO;
ir   = iresponse(Phi,Sigma,hor,Omeg);


end
