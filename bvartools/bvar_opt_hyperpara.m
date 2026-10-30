function minus_log_dnsty = bvar_opt_hyperpara(hyperpara,y,lags,options)
%BVAR_OPT_HYPERPARA Objective for BVAR_MAX_HYPER: minus the log marginal likelihood.
%
%   minus_log_dnsty = BVAR_OPT_HYPERPARA(hyperpara,y,lags)
%   minus_log_dnsty = BVAR_OPT_HYPERPARA(hyperpara,y,lags,options)
%
%   Inputs
%     hyperpara  log of the hyperparameters being estimated
%     y, lags    as in BVAR_ML
%     options    as in BVAR_ML, plus .index_est, .index_fixed and
%                .hyperpara_fidex (log values of the fixed ones), set by
%                BVAR_MAX_HYPER
%
%   Output
%     minus_log_dnsty   scalar
%
%   This is the default objective of BVAR_MAX_HYPER.
%
%   See also BVAR_MAX_HYPER, BVAR_ML.

%   Filippo Ferroni, 6/1/2017.


HYPERPARA_LB = 1e-5;
HYPERPARA_UB = 1e4;



if nargin <4
    hyperpara = exp(hyperpara);


    if any(~isfinite(hyperpara)) || any(~isreal(hyperpara)) || ...
            any(hyperpara < HYPERPARA_LB) || any(hyperpara > HYPERPARA_UB)
        minus_log_dnsty = 1e10;
        return
    end

    log_dnsty = bvar_ml(hyperpara,y,lags);
else
    hyperpara(options.index_est)   = exp(hyperpara);
    if isempty(options.index_fixed) == 0
        hyperpara(options.index_fixed) = exp(options.hyperpara_fidex);
    end


    if any(~isfinite(hyperpara)) || any(~isreal(hyperpara)) || ...
            any(hyperpara < HYPERPARA_LB) || any(hyperpara > HYPERPARA_UB)
        minus_log_dnsty = 1e10;
        return
    end



    if isfield(options,'pandemic_priors') == 1 
        log_dnsty = bvar_ml_pandemic(hyperpara,y,lags,options,options.pandemic_priors) ;
    else
        log_dnsty                      = bvar_ml(hyperpara,y,lags,options);
    end



end

minus_log_dnsty = -log_dnsty;