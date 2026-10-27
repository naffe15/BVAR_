function out = fast_ols_(y, X)
%FAST_OLS_ OLS by singular value decomposition.
%
%   out = FAST_OLS_(y,X)
%
%   Inputs
%     y   T x ny   dependent variables
%     X   T x nk   regressors
%
%   Output
%     out.B     nk x ny   coefficients
%     out.u     T x ny    residuals
%     out.xxi   nk x nk   inv(X'X), computed from the SVD
%     out.y, out.X        the inputs
%
%   See also RFVAR3, YXB_.

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
