function y = ivech(x)
%IVECH Symmetric matrix from its VECH.
%
%   y = IVECH(x)
%
%   Input    x   n(n+1)/2 x 1   output of VECH
%   Output   y   n x n          symmetric matrix
%
%   See also VECH.

[K,k] = size(x);
if  k ~= 1
    error('x must be a column vector');
end

N = -1/2 + sqrt(1/4+2*K);
if floor(N) ~= N
    error('something went wrong');
end

y = nan(N,N);
a = 1; 
for jj = 1 : N
    b = a + N - jj; 
    span = a : b;
    y(jj,jj:N) = x(span,1)'; % row
    y(jj:N,jj) = x(span,1); % column
    a = b + 1; 
end

