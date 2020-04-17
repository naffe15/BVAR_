function log_dnsty = blp_ml(hyperpara,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'blp_ml' computes the marginal likelihood for the NMIW

% Inputs:
% - hyperpara, shrinkage hyperpara over which maximize the marginal
% likelihood

% Output: marginal data density

% Filippo Ferroni, 3/21/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************
% SETTINGS
%********************************************************
ny = size(G,2);
%********************************************************
% Conjugate Prior: N-IW
%********************************************************
% constructing the prior mean

prior.tau(hh)     = hyperpara; 
[posterior,prior] = p2p(hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant);


prior_int     = matrictint(prior.S, prior.df, prior.XXi);
posterior_int = matrictint(posterior.S, posterior.df, posterior.XXi);
lik_nobs      = posterior.df - prior.df;
log_dnsty     = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);

