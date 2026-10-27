function [ir,Omeg] = iresponse_sign_hmoments_eig(errors,Phi,Sigma,hor,moment)
%IRESPONSE_SIGN_HMOMENTS_EIG Responses identified from the higher moments of the residuals.
%
%   [ir,Omeg] = IRESPONSE_SIGN_HMOMENTS_EIG(errors,Phi,Sigma,hor)
%   [ir,Omeg] = IRESPONSE_SIGN_HMOMENTS_EIG(errors,Phi,Sigma,hor,moment)
%
%   The rotation is taken from the singular value decomposition of the
%   fourth-order cumulant (moment = 4) or of the third moments (moment = 3)
%   of the residuals. Used by BVAR_ with options.hmoments_eig.
%
%   Inputs
%     errors   T x ny    reduced-form residuals
%     Phi      (ny*lags+nx) x ny   VAR coefficients (see IRESPONSE)
%     Sigma    ny x ny   covariance of the reduced-form residuals
%     hor      scalar    number of horizons
%     moment   3 or 4    moments used (default 4)
%
%   Outputs
%     ir       ny x hor x ny   responses
%     Omeg     ny x ny         rotation of chol(Sigma,'lower')
%
%   See also IRESPONSE_SIGN_HMOMENTS, THIRDMOM, FOURTHMOM.

%   Filippo Ferroni.

[~,n]   = size(Sigma);
if nargin < 5
    moment  = 4;    
end

A       = chol(Sigma,'lower');
%  varols.u * BVAR.iSigChol';
% u     = errors * A';   % old: does not standardize the residuals
u     = errors / A';     % standardized residuals, u = errors * inv(A)', cov(u) = I

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
