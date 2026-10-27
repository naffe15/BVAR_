function [x,info]=lyapunov_symm(a,b)
%LYAPUNOV_SYMM Solves the symmetric Lyapunov equation x - a*x*a' = b.
%
%   [x,info] = LYAPUNOV_SYMM(a,b)
%
%   Used for the unconditional covariance of a stationary state vector.
%
%   Inputs
%     a     n x n   transition matrix (stable)
%     b     n x n   symmetric matrix
%
%   Outputs
%     x     n x n   symmetric solution
%     info  0
%
%   See also KFILTERNAN, P2P.

%   Copyright (C) 2006-2017 Dynare Team. Dynare is free software released
%   under the GNU General Public License, version 3 or later.

  info = 0;
  n = size(b,1);
  if n == 1
    x=b/(1-a*a);
    return
  end
  x=zeros(n,n);
  [u,t]=schur(a);
  b=u'*b*u;
  for i=n:-1:2
    if t(i,i-1) == 0
      if i == n
	c = zeros(n,1);
      else
	c = t(1:i,:)*(x(:,i+1:end)*t(i,i+1:end)')+...
	    t(i,i)*t(1:i,i+1:end)*x(i+1:end,i);
      end
      q = eye(i)-t(1:i,1:i)*t(i,i);
      x(1:i,i) = q\(b(1:i,i)+c);
      x(i,1:i-1) = x(1:i-1,i)';
    else
      if i == n
	c = zeros(n,1);
	c1 = zeros(n,1);
      else
	c = t(1:i,:)*(x(:,i+1:end)*t(i,i+1:end)')+...
	    t(i,i)*t(1:i,i+1:end)*x(i+1:end,i)+...
	    t(i,i-1)*t(1:i,i+1:end)*x(i+1:end,i-1);
	c1 = t(1:i,:)*(x(:,i+1:end)*t(i-1,i+1:end)')+...
	     t(i-1,i-1)*t(1:i,i+1:end)*x(i+1:end,i-1)+...
	     t(i-1,i)*t(1:i,i+1:end)*x(i+1:end,i);
      end
      q = [eye(i)-t(1:i,1:i)*t(i,i) -t(1:i,1:i)*t(i,i-1);...
	   -t(1:i,1:i)*t(i-1,i) eye(i)-t(1:i,1:i)*t(i-1,i-1)];
      z =  q\[b(1:i,i)+c;b(1:i,i-1)+c1];
      x(1:i,i) = z(1:i);
      x(1:i,i-1) = z(i+1:end);
      x(i,1:i-1)=x(1:i-1,i)';
      x(i-1,1:i-2)=x(1:i-2,i-1)';
      i = i - 1;
    end
  end
  if i == 2
    c = t(1,:)*(x(:,2:end)*t(1,2:end)')+t(1,1)*t(1,2:end)*x(2:end,1);
    x(1,1)=(b(1,1)+c)/(1-t(1,1)*t(1,1));
  end
  x=u*x*u';