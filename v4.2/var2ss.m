function [A,B,C,const,Sigma,lags,index_var]=var2ss(Phi,Sigma,index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% var2ss coverts the VAR into a state space system of this form
% x(t) = A x(t-1) + B Sigma' u(t) ~ N(0,I)
% y(t) = C*(cons + x(t-1))
% where A is the companion form of the lag struture

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    %% all stocks
    index = zeros(size(Sigma,1),1);
end

Sigma = chol(Sigma);

N           = size(Sigma,1);
[m , n]     = size(Phi);
lags        = (m-1)/n;

% companion form
A       = [Phi(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
B       = eye(N * lags, N);
C       = eye(N, N * lags);

IminusAlags = eye(N);
for ell = 1 : lags
    IminusAlags = IminusAlags - Phi(N * (ell -1) + 1 : N * ell, :)'; 
end
iIminusAlags  = inv(IminusAlags);

const    = [iIminusAlags*Phi(end, :)'; zeros(N*(lags-1), 1)];

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

