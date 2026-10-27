function [A,B,C,const,Sigma,index_var]=uc2ss(Phi,Sigma)
%UC2SS State-space form of the unobserved-components model of UC2_.
%
%   [A,B,C,const,Sigma,index_var] = UC2SS(Phi,Sigma)
%
%   Model:
%     y(t) = a(t) + b(t)
%     a(t) = a1 a(t-1) + ... + ap a(t-p) + ea(t)   (cycle)
%     b(t) = c(t-1) + b(t-1) + eb(t)               (trend)
%     c(t) = c(t-1) + ec(t)                        (trend growth)
%   written as  x(t) = A x(t-1) + B Sigma' u(t),  y(t) = C (const + x(t)).
%
%   Inputs
%     Phi     p x 1   autoregressive coefficients of the cycle
%     Sigma   variances of (ea, eb, ec)
%
%   Outputs
%     A, B, C, const   state-space matrices (const = 0)
%     Sigma            standard deviations of the shocks
%     index_var        state of the observable
%
%   See also UC2_, KFILTERNAN, VAR2SS.

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