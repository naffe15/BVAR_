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
constant_ = 0;
if size(alpha,1) > N * lags
    constant_ = 1;
end
exogenous_ = 0;
nexogenous = 0;
if size(alpha,1) > N * lags + 1
    exogenous_ = 1;
    exogenous  = bvar_.XX(:,N*lags + 2 : end);
    nexogenous = size(exogenous,2);
end

% Retrieve the initial condition
yo     = bvar_.XX(1,1:N*lags)'; % remove the constant, timetrend or exogenous variables from X

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

Tu     = length(u);
A      = chol(Sigma,'lower');
ierror = zeros(Tu,N); 

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
C       = 0;    
% if constant
if constant_
    C = [alpha(N * lags +1, :)'; zeros(N*(lags-1),1)];
end
% if exogenous
if exogenous_
    Gamma = [alpha(N * lags + 2 : end, :)'; zeros(N*(lags-1),nexogenous)];
end

    
%ystar   = zeros(length(u),N*lags,N);

% Deterministic Part
B_ = zeros(Tu,N);
for t = 1 : Tu
    Aa =  0 ;
    for tau = 1 : t
        Aa = Aa + F^(tau-1)*C;
    end
%     B(t,:) = (Aa + F^(t) * yo)';
    B  = Aa + F^(t) * yo;
    B_(t,:) = B(1 : N,1)';
%     initial condition
%    Bb(t,:) = (F^(t) * yo)';
end

Kappa = G * A * Omega ;
% Stochastic part
E_ = zeros(Tu,N,N);
for shock = 1 : N
    Ind              = zeros(N); 
    Ind(shock,shock) = 1;
    for t = 1 : Tu
        D_  = 0;
        for tau = 1 : t
            D_     = D_ + F^(t-tau) *  Kappa * Ind * ierror(tau,:)';
        end
        E_(t,:,shock) = D_(1:N,1)';
    end
end

% Exogenous Part (if any)
Q_ = 0;
if exogenous_
    Q_ = zeros(Tu,N,nexogenous);
    for exo = 1 : nexogenous
        for t = 1 : length(u)
            M_  = 0;
            for tau = 1 : t
                M_     = M_ + F^(t-tau) *  Gamma(:,exo) * exogenous(tau,exo)';
            end
            Q_(t,:,exo) = M_(1:N,1)';
        end
    end
end


% check
W       = B_ + sum(E_,3) + sum(Q_,3);
yhat    = W(:, 1 : N);
tmp = max(max(abs(data(lags+1:end,:) - yhat)));
if  tmp> tol
    histdec = [];
    warning(['Maximium Discrepancy ' num2str(tmp)])
    return;
end


histdec             = E_; % stochastic part
if exogenous_
    histdec(:,:, N+1 : N+nexogenous)    = Q_; % exogenous part
end
if constant_ 
    histdec(:,:, N+nexogenous+1 )    = B_; % deterministic part
end


end

