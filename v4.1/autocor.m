function A = autocor(X,m)
%AUTOCOR     Computes cross-correlations and autocorrelations.
%            A = autocor(X,m) takes an T x N matrix X and constructs
%
%                           [ B[0] ]
%                           [ B[1] ]
%                    A    = [ B[2] ]
%                           [ B[3] ]
%                           [  :   ]
%                           [ B[m] ]
%
%            where the (i,j)th element of B[k] is corr(X[i,t],X[j,t-k]),
%            i,j in [1,...,N].  A has dimension N*(m+1) x N.
%
 
%            Ellen R. McGrattan, 2-1-99
%
[T,N]=size(X);
if T-m<=0;
  str='The number of obs.in the time series must be larger than m';
  error(str);
end;
b = zeros(N);
A = zeros(N*(m+1),N);
for k=1:m+1;
  for i=1:N;
    for j=1:N;
      x      = corrcoef([X(k:T,i),X(1:T-k+1,j)]);
      b(i,j) = x(2,1);
    end
  end
  A((k-1)*N+1:k*N,:) = b;
end
 
 
