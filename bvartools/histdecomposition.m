function [histdec,ierror] = histdecomposition(error,alpha,Sigma,data,Omega,yo,tol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'histdecomposition' computes decomposition of the observable variables in
% terms of structural VAR shocks and initial condition

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin < 7
    tol =  10^(-7);    
end

N     = size(Sigma,1);
[m , k]     = size(alpha);
lags        = (m-1)/k;

if nargin < 5    
    Omega = eye(N);
    yo    = zeros(N*lags,1);
end
if nargin < 6   
    yo    = zeros(N*lags,1);
end
% 

A      = chol(Sigma,'lower');
ierror = zeros(length(error),N); 

% (1) e = A * Omega * eta 
% e = n * 1         %reduced 
% eta = n * 1       % structural
% A = n * n 
% Omega = n * n
% e' =  eta' * A' * Omega';
% (2) E  =  ETA  * A' * Omega';
% where
% H = T * n
% E = T * n



ierror = error * inv( Omega' * A'); %#ok<MINV>


% companion
F       = [alpha(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
C       = [alpha(end, :)'; zeros(N*(lags-1),1)];
ystar   = zeros(length(error),N*lags,N);

% Deterministic Part
for t = 1 : length(error)
    Aa =  0 ;
    for tau = 1 : t
        Aa = Aa + F^(tau-1)*C;
    end
%     B(t,:) = (Aa + F^(t) * yo)';
    B(t,:) = (Aa + F^(t) * yo)';
%     initial condition
    Bb(t,:) = (F^(t) * yo)';
end

Kappa = G * A * Omega ;
% Stochastic part
for shock = 1 : N
    Ind              = zeros(N); 
    Ind(shock,shock) = 1;
    for t = 1 : length(error)
        D  = 0;
        for tau = 1 : t
            D     = D + F^(t-tau) *  Kappa * Ind * ierror(tau,:)';
        end
        E(t,:,shock) = D';
    end
end

% check
W       = B + sum(E,3);
yhat    = W(:, 1 : N);
tmp = max(max(abs(data(lags+1:end,:) - yhat)));
if  tmp> tol,
    histdec = [];
    warning(['Maximium Discrepancy ' num2str(tmp)])
    return;
end

histdec             = E(:,1:N,:);
histdec(:,:,N+1)    = B(:,1:N);

end

