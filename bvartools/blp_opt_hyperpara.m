function minus_log_dnsty = blp_opt_hyperpara(hyperpara,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant)
%BLP_OPT_HYPERPARA Objective for the shrinkage of a Bayesian local projection.
%
%   minus_log_dnsty = BLP_OPT_HYPERPARA(hyperpara,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant)
%
%   Passed by name to the optimiser in DIRECTMETHODS.
%
%   Inputs
%     hyperpara   log of the shrinkage
%     others      as in P2P
%
%   Output
%     minus_log_dnsty   minus the log marginal likelihood
%
%   See also BLP_ML, DIRECTMETHODS.

%   Filippo Ferroni, 3/21/2020.

hyperpara       = exp(hyperpara);
log_dnsty       = blp_ml(hyperpara,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant);
minus_log_dnsty = -log_dnsty;