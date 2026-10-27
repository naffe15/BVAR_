function y = vech(x)
%VECH Stacks the lower triangle of a square matrix, column by column.
%
%   y = VECH(x)
%
%   Input    x   n x n
%   Output   y   n(n+1)/2 x 1
%
%   See also IVECH, DUPLICATIONMATRIX.

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
