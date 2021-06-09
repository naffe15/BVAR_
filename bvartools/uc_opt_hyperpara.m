function minus_log_dnsty = uc_opt_hyperpara(hyperpara,y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% computes the likelihood of the US system
% a(t)  = a1 a(t-1) + ... + ap a(t-p) + ea(t )    [transition 1]
% b(t)  = c(t-1)    + b(t-1) + eb(t);             [transition 2]
% c(t)  = c(t-1)             + ec(t);             [transition 3]
% y(t)  = a(t)      + b(t);                       [transistion 4]

% to a state space of the form
% x(t) = A x(t-1) + B Sigma' u(t) ~ N(0,I)
% y(t) = C*(cons + x(t-1))
% where A is the companion form of the lag struture

% Filippo Ferroni, 6/1/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% options.state_space_model = 2;
% 
% options.only_logL    = 1; 
% options.initialCond  = 2; 
% options.aZero        = zeros(lags+2+1,1);
% options.aZero(end-2) = y(1);
% options.pZero        = 10*eye();

hyperpara(options.index_est)   = exp(hyperpara);
if isempty(options.index_fixed) == 0
    hyperpara(options.index_fixed) = exp(options.hyperpara_fidex);
end

phi   = (hyperpara(1:lags));
sigma = diag((hyperpara(1+lags:3+lags)));

log_dnsty       = kfilternan(phi,sigma,y,options);
minus_log_dnsty = -log_dnsty;