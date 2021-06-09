function [histdec,ierror] = histdecomp(bvar_,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'histdecomp' computes decomposition of the observable variables in
% terms of structural VAR shocks and initial condition

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tolerance
tol =  10^(-7);
% Retrieve the initial condition
yo     = bvar_.XX(1,1:end-1)'; % remove the constant from X
% use posterior means
u      = mean(bvar_.e_draws,3);       % reduced form errors
alpha  = mean(bvar_.Phi_draws,3);     % AR
Sigma  = mean(bvar_.Sigma_draws,3);   % Sigma
% retrieve VAR settings
N           = bvar_.N;
lags        = bvar_.lags;
% no rotation
Omega = eye(N);
% data
data = bvar_.data;


if nargin > 1
    if isfield(options,'tol')==1
        tol = options.tol;
    end
    if isfield(options,'yo')==1 && options.yo == 0 %no initial condition
        yo    = zeros(N*lags,1);
    end
    if isfield(options,'Omega')==1 % use a specific rotation
         Omega = options.Omega;
    end
    if isfield(options,'draw')==1 % use a specific draw for the reduced form param
         draw = options.draw;
         u      = bvar_.e_draws(:,:,draw);            % reduced form errors
         alpha  = bvar_.Phi_draws(:,:,draw);          % AR
         Sigma  = bvar_.Sigma_draws(:,:,draw);        % Sigma
    end
    if isfield(options,'median')==1 && options.median == 1 % use median instead of mean
        u      = median(bvar_.e_draws,3);       % reduced form errors
        alpha  = median(bvar_.Phi_draws,3);     % AR
        Sigma  = median(bvar_.Sigma_draws,3);   % Sigma
    end    
end

A      = chol(Sigma,'lower');
ierror = zeros(length(u),N); 

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

% structural innovations
ierror = u * inv( Omega' * A'); %#ok<MINV>

% companion form
F       = [alpha(1 : N * lags, :)'; eye(N*(lags-1), N*lags)];
G       = eye(N * lags, N);
C       = [alpha(end, :)'; zeros(N*(lags-1),1)];
ystar   = zeros(length(u),N*lags,N);

% Deterministic Part
for t = 1 : length(u)
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
    for t = 1 : length(u)
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

