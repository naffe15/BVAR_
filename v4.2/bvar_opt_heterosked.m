function minus_log_dnsty = bvar_opt_heterosked(hyperpara,y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bvar_opt_heterosked_' computes the marginal likelihood over the
% hyperparameters of for the Minnesota prior and the rescaling parameters, 
% see Lenza and Primiceri (2020)   

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