function minus_log_dnsty = bvar_opt_heterosked(hyperpara,y,lags,options)
%BVAR_OPT_HETEROSKED Objective for the Minnesota and COVID rescaling hyperparameters.
%
%   minus_log_dnsty = BVAR_OPT_HETEROSKED(hyperpara,y,lags,options)
%
%   Pass it to BVAR_MAX_HYPER with options.objective_function =
%   'bvar_opt_heterosked'.
%
%   Inputs
%     hyperpara  log of the estimated hyperparameters: the five Minnesota
%                ones, then one scale for each period from options.tstar on
%     y, lags    as in BVAR_ML
%     options    as in BVAR_OPT_HYPERPARA, plus .tstar, the first period
%                with a rescaled residual variance
%
%   Output
%     minus_log_dnsty   scalar, including the Jacobian of the rescaling
%
%   See also BVAR_MAX_HYPER, BVAR_ML.
%
%   Reference: Lenza, M. and G. E. Primiceri (2022), "How to Estimate a
%   Vector Autoregression after March 2020", Journal of Applied
%   Econometrics 37(4), 688-699.

%   Filippo Ferroni, 6/1/2017.

hyperpara(options.index_est)   = exp(hyperpara);
if isempty(options.index_fixed) == 0
    hyperpara(options.index_fixed) = exp(options.hyperpara_fidex);
end

esse   = ones(size(y,1),1);
nhpara = length(hyperpara);
% hyper para for covid
% s0 (at tstar), s1, s2 and rho
jj = 0;
for nh = 6 : nhpara % for nh = 6 : nhpara-1
    esse(options.tstar + nh-6,:)   = hyperpara(nh);
    jj = jj+1;
end
% jj0 = jj;
% for nj = options.tstar + 3 : size(y,1)
%     jj = jj+1;
%     esse(nj ,:)   = 1 + (1-hyperpara(end-1))*hyperpara(end)^(jj-jj0);
% end
esse(1:lags) = []; 
options.heterosked_weights = esse;

log_dnsty       = bvar_ml(hyperpara,y,lags,options);
log_dnsty       = log_dnsty - size(y,2) * sum(log(esse(lags+1:end,1)));
minus_log_dnsty = -log_dnsty;