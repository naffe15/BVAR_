% function [log_dnsty] = blp_ml(shrinkage,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant)
function [log_dnsty] = blp_ml(shrinkage,hh,prior,olsreg_,F,G,Fo,positions_nylags,position_constant)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'blp_ml' computes the marginal likelihood for the NMIW LP

% Inputs:
% - hyperpara, shrinkage hyperpara over which maximize the marginal
% likelihood

% Output: marginal data density

% Filippo Ferroni, 3/21/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************
% SETTINGS
%********************************************************
ny   = size(G,2);
lags = size(G,2)/ny;
%********************************************************
% Conjugate Prior: MN-IW for LP
%********************************************************
[posterior1,prior1] = p2p(hh,shrinkage,prior,olsreg_,F,G,Fo,positions_nylags,position_constant);

%*******************************************************************
%* Compute the log marginal data density for the VAR model with MNIW
%*******************************************************************
var.y = olsreg_.Y;
var.X = olsreg_.X(:,[positions_nylags,position_constant]);
prior1.Sigma.df     = prior1.df;
prior1.Sigma.scale  = prior1.S;
prior1.Phi.mean     = prior1.BetaMean;
prior1.Phi.cov      = prior1.BetaVar;

log_dnsty      = mniw_log_dnsty(prior1,posterior1,var);

