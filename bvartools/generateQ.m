function Omega = generateQ(m)

% generate an orthonormal matrix.

G     = randn(m);
[Q,R] = qr(G);
% normalize to positive entry in the diagonal
In    = diag(sign(diag(R)));
Omega = Q  * In;