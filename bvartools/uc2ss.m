function [A,B,C,const,Sigma,index_var]=uc2ss(Phi,Sigma)

% from the UC of the form 
% y(t)  = a(t)      + b(t);                       [transition 0]
% a(t)  = a1 a(t-1) + ... + ap a(t-p) + ea(t )    [transition 1]
% lags companion
% b(t)  = c(t-1)    + b(t-1) + eb(t);             [transition 2]
% c(t)  = c(t-1)             + ec(t);             [transition 3]

% to a state space of the form
% x(t) = A x(t-1) + B Sigma' u(t) ~ N(0,I)
% y(t) = C*(cons + x(t-1))
% where A is the companion form of the lag struture

if size(Phi,2) > size(Phi,1)
    Phi = Phi'; % autoregressive parameters for the cycle
end

N    = 1;
lags = length(Phi);
% companion form cycle
A0     = [Phi'; eye(N*(lags-1), N*lags)];
% companion form trend
A1     = [1 1; 0  1];
% 
A = [A0 zeros(size(A0,1),2);
    zeros(2,size(A0,2)) A1];
C           = zeros(1,size(A,1));
C(1,1)      = 1;
C(1,end-1)  = 1;
A = [A; C];
A = [A, zeros(size(A,1),1)];

index_var = size(A,1);

B     = zeros(index_var,3);
B(1,1)         = 1; % cycle
B(end-1,end-1) = 1; % trend1
B(end,end)     = 1; % trend2

C           = zeros(1,index_var);
C(1,end)    = 1;

const = zeros(index_var,1);
Sigma = std(Sigma);