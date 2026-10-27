function [A,B,C,const,Sigma,lags,index_var]=var2ss(Phi,Sigma,index)
%VAR2SS State-space form of a VAR, with temporal aggregation.
%
%   [A,B,C,const,Sigma,lags,index_var] = VAR2SS(Phi,Sigma)
%   [A,B,C,const,Sigma,lags,index_var] = VAR2SS(Phi,Sigma,index)
%
%   x(t) = A x(t-1) + B Sigma' u(t),   u(t) ~ N(0,I)
%   y(t) = C (const + x(t))
%   where A is the companion matrix, possibly with extra states for the
%   aggregated variables.
%
%   Inputs
%     Phi     (ny*lags+nx) x ny   VAR coefficients (nx = 0 or 1, constant)
%     Sigma   ny x ny   residual covariance
%     index   ny x 1    how each variable is observed (default 0):
%               0  stock, observed as is
%               2  flow or deflator: average of the last 3 periods
%                  (monthly to quarterly); needs lags >= 2
%               4  average of the last 4 periods; needs lags >= 3
%
%   Outputs
%     A, B, C     state-space matrices
%     const       steady state of the states
%     Sigma       ny x ny  upper Cholesky factor of the input Sigma
%     lags        number of lags
%     index_var   1 x ny   state that each observable loads on
%
%   See also KFILTERNAN, UC2SS, CFORECASTS.

%   Filippo Ferroni, 6/1/2015.

if nargin < 3
    %% all stocks
    index = zeros(size(Sigma,1),1);
end

Sigma = chol(Sigma);

N           = size(Sigma,1);
[m , n]     = size(Phi);
if rem(m, n)==0 % no constant
    lags = m/n;
    nx   = 0;
else % with constant
    nx  =  1;
    lags = floor((m-1)/n);
end

% companion form
A       = [Phi(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
B       = eye(N * lags, N);
C       = eye(N, N * lags);

IminusAlags = eye(N);
for ell = 1 : lags
    IminusAlags = IminusAlags - Phi(N * (ell -1) + 1 : N * ell, :)'; 
end
%iIminusAlags  = inv(IminusAlags);
iIminusAlags  = IminusAlags\eye(N);

if nx==1
    const    = [iIminusAlags*Phi(end, :)'; zeros(N*(lags-1), 1)];
else
    const    = [iIminusAlags*zeros(N,1); zeros(N*(lags-1), 1)];
end

% index =0 % stock:      xq(t) = xm(t)
% index =1 % TBA
% index =2 % deflator/real flow:   xq(t) = 1/3( xm(t) +  xm(t-1) +  xm(t-2))

index_var = 1 : N;

for vv = 1 : N    
    if index(vv) == 1 
%         C(vv,vv : N : N * 3) = 1;
    elseif index(vv) == 2 % 
        if lags < 2
            error('When you specify the flow/deflation aggregation, you need at least 2 lags')
        end
        A                     = [A zeros(size(A,1),1)];
        A                     = [A; zeros(1,size(A,2))];
        A(end,vv : N : N * 2) = 1/3;
        Ao                    = eye(length(A));
        Ao(end,vv)            = -1/3;
        iAo                   = inv(Ao);
        A                     = iAo * A;
        
        B = [B; zeros(1,size(B,2))];
        B = iAo * B;
        
        C  = [C zeros(N,1)];
        C(vv,vv)  = 0;
        C(vv,end) = 1;
        
        const(size(A,1)) = const(vv);
             
        index_var(vv) = size(A,1);
        
    elseif index(vv) == 4 % 
        if lags < 3
            error('When you specify the flow/deflation aggregation (weekly-monthly or quarterly-annual), you need at least 3 lags')
        end
        A                     = [A zeros(size(A,1),1)];
        A                     = [A; zeros(1,size(A,2))];
        A(end,vv : N : N * 3) = 1/4;
        Ao                    = eye(length(A));
        Ao(end,vv)            = -1/4;
        iAo                   = inv(Ao);
        A                     = iAo * A;
        
        B = [B; zeros(1,size(B,2))];
        B = iAo * B;
        
        C  = [C zeros(N,1)];
        C(vv,vv)  = 0;
        C(vv,end) = 1;
        
        const(size(A,1)) = const(vv);
             
        index_var(vv) = size(A,1);
        
    end

end

