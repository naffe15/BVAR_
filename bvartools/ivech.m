function y = ivech(x)
% inverse of vech (vectorize the lower triangular part of a matrix x)

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

