function out = fast_ols_(y, X)
% FAST_OLS_  SVD-based OLS: returns coefficients, residuals, and (X'X)^{-1}.
%
%   out = fast_ols_(y, X)
%
%   Inputs
%     y  -- (T x ny) matrix of dependent variables
%     X  -- (T x nk) matrix of regressors
%
%   Output fields
%     out.B    -- (nk x ny) OLS coefficient matrix
%     out.u    -- (T  x ny) residuals
%     out.xxi  -- (nk x nk) (X'X)^{-1} computed via SVD
%     out.y    -- y  (passed through)
%     out.X    -- X  (passed through)

[vl, d_, vr] = svd(X, 0);
di  = 1./diag(d_);
nk  = size(X, 2);
B   = (vr .* repmat(di', nk, 1)) * vl' * y;
u   = y - X * B;
xxi = vr .* repmat(di', nk, 1);
xxi = xxi * xxi';

out.B   = B;
out.u   = u;
out.xxi = xxi;
out.y   = y;
out.X   = X;
end
