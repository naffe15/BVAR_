function minus_log_dnsty = uc_opt_hyperpara(hyperpara,y,lags,options)
%UC_OPT_HYPERPARA Objective for UC2_: minus the log likelihood of the UC model.
%
%   minus_log_dnsty = UC_OPT_HYPERPARA(hyperpara,y,lags,options)
%
%   Inputs
%     hyperpara  estimated parameters as seen by the optimiser: the lags AR
%                coefficients of the cycle in levels, then the three shock
%                variances in logs
%     y          T x 1    series
%     lags       scalar   order of the cycle
%     options    struct   .index_est, .index_fixed, .hyperpara_fidex (fixed
%                         values, same scale as hyperpara), plus the
%                         KFILTERNAN options
%
%   Output
%     minus_log_dnsty   scalar
%
%   See also UC2_, UC2SS, KFILTERNAN.

%   Filippo Ferroni, 6/1/2020.

% options.state_space_model = 2;
% 
% options.only_logL    = 1; 
% options.initialCond  = 2; 
% options.aZero        = zeros(lags+2+1,1);
% options.aZero(end-2) = y(1);
% options.pZero        = 10*eye();

% hyperpara(options.index_est)   = exp(hyperpara);   % old: everything in logs, so
% if isempty(options.index_fixed) == 0                 % the AR coefficients were
%     hyperpara(options.index_fixed) = exp(options.hyperpara_fidex);   % forced positive
% end
xopt      = hyperpara;
hyperpara = zeros(lags+3,1);
hyperpara(options.index_est) = xopt;
if isempty(options.index_fixed) == 0
    hyperpara(options.index_fixed) = options.hyperpara_fidex;
end
islog            = (1:lags+3)' > lags;    % AR coefficients (1:lags) in levels, variances in logs
hyperpara(islog) = exp(hyperpara(islog));

phi   = (hyperpara(1:lags));
sigma = diag((hyperpara(1+lags:3+lags)));

log_dnsty       = kfilternan(phi,sigma,y,options);
minus_log_dnsty = -log_dnsty;