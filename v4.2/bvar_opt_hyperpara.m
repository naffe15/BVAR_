function minus_log_dnsty = bvar_opt_hyperpara(hyperpara,y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bvar_opt_hyperpara' computes the marginal likelihood over the
% hyperparameters of for the Minnesota prior  
% Bridge function between 'bvar_ml' and the minimization routine
% 'bvar_max_hyper' 

% Inputs:
% - hyperpara, Minnesota hyperpara over which maximize the marginal
% likelihood
% - y, data columns variables
% - lags, lag order of the VAR
% - options, see below for details

% Output: marginal data density 

% Filippo Ferroni, 6/1/2017
% Revised, 3/21/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin <4
    hyperpara = exp(hyperpara);
    log_dnsty = bvar_ml(hyperpara,y,lags);
else
    hyperpara(options.index_est)   = exp(hyperpara);
    if isempty(options.index_fixed) == 0
        hyperpara(options.index_fixed) = exp(options.hyperpara_fidex);
    end
    log_dnsty                      = bvar_ml(hyperpara,y,lags,options);
end

minus_log_dnsty = -log_dnsty;