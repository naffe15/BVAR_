function minus_log_dnsty = blp_opt_hyperpara(hyperpara,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'blp_opt_hyperpara' computes the marginal likelihood over the
% hyperparameters of for the LP  

% Inputs:
% - hyperpara, prior shrinkage
% likelihood

% Output: marginal data density 

% Filippo Ferroni, 3/21/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hyperpara       = exp(hyperpara);
log_dnsty       = blp_ml(hyperpara,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant);
minus_log_dnsty = -log_dnsty;