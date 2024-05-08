function Omega = generateQ(m,N)

% generate an orthonormal matrix.

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