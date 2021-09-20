function [Sig,Kap] = GetSigKap(A,C,Sigma1,Sigma2)
%GetSigKap computes the steady state Kalmna gain
% x(t) = A x(t-1) + e1(t) ~ (0,Sigma1)
% y(t) = C x(t-1) + e2(t) ~ (0,Sigma2)

ny = size(C,1);
ns = size(A,1);

% starting values
Sig0 = zeros(ns);
Kap0 = zeros(ns,ny);

diff = 1;
% set the tolerance level for the convergence
tol = 1E-6;
JJ  = 1E6; jj = 0; 

while diff > tol
    %Kap = (A*Sig0*C')*inv(C*Sig0*C' + Sigma2);
    Kap  = (A*Sig0*C')/(C*Sig0*C' + Sigma2);
    Sig  = A*Sig0*A' + Sigma1 - Kap0 * C*Sig0*A';
    %diff = max([max( abs(Kap-Kap0))  abs( vec(Sig-Sig0)) ]);
    diff = max([max(max( abs(Kap-Kap0) )) max(max( abs(Sig-Sig0) ))])  ;
    Sig0 = Sig; Kap0 = Kap;
    jj   = jj + 1;
    if jj == JJ
        error('Could not compute the steady state Kalman Gain')
    end    
end

end