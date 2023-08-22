function y = vech(x)
% vectorize the lower triangular part of a matrix x
if size(x,1) ~= size(x,2)
    error('x must be a square matrix');
end
N = size(x,2);
y = [];%nan(N*(N+1)/2);
for jj = 1 : N
    y = [y; x(jj:end,jj)];
end


% y = reshape(tril(x),size(x,1)*size(x,2),1);
% y(y ==0) = [];
% if issymmetric(x) == 0
%     error('x must be a symmetrix matrix');
% end
