function Omega = generateQ(m,N)
%GENERATEQ Random orthonormal matrices, uniformly distributed (Haar).
%
%   Omega = GENERATEQ(m)
%   Omega = GENERATEQ(m,N)
%
%   Inputs
%     m       scalar      dimension of the matrix
%     N       scalar      number of matrices (default 1)
%
%   Output
%     Omega   m x m x N   orthonormal matrices from the QR decomposition of
%                         a standard normal matrix, with the diagonal of R
%                         normalised to be positive
%
%   See also IRESPONSE_SIGN, IRESPONSE_SET, MAX_FEVD.
%
%   Reference: Rubio-Ramirez, J. F., D. F. Waggoner and T. Zha (2010),
%   Review of Economic Studies 77(2), 665-696.

if nargin<2
    N = 1;
    Omega = nan(m);
end

for n = 1 : N
    G     = randn(m);
    [Q,R] = qr(G);
    % normalize to positive entry in the diagonal
    In    = diag(sign(diag(R)));
    Omega(:,:,n) = Q  * In;
end