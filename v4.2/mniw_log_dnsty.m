function  log_dnsty = mniw_log_dnsty(prior,posterior,var)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'logmlike' computes the marginal likelihood for the NM-IW 
% Inputs:
% Output: marginal data density
% Filippo Ferroni, 3/21/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nobs,ny]  = size(var.y);
nk         = size(var.X,2);

Fv    = chol(prior.Phi.cov)';
Fo    = chol(inv(prior.Sigma.scale))';
iV    = inv(prior.Phi.cov);

var.e = var.y - var.X*posterior.PhiHat;

log_dnsty = - nobs*ny/2 * log(pi);
log_dnsty = log_dnsty - nobs / 2 *log(det(prior.Sigma.scale));
log_dnsty = log_dnsty - ny/2 * log(det( eye(nk) + Fv'*var.X'*var.X*Fv ));
% % Giannone, Lenza Primiceri (2015) Appendix  A13-A14
FF    = eye(ny) ...
    + Fo'* ((posterior.PhiHat - prior.Phi.mean)'* iV * (posterior.PhiHat - prior.Phi.mean) ...
    + var.e'* var.e) * Fo;
log_dnsty = log_dnsty - (nobs + prior.Sigma.df)/2 * log(det(  FF ));
log_dnsty = log_dnsty + ggammaln(ny,(nobs + prior.Sigma.df)/2) ;
log_dnsty = log_dnsty - ggammaln(ny,prior.Sigma.df/2);


function lgg = ggammaln(m, df)
if df <= (m-1)
    error('too few df in ggammaln; increase the number of observations or reduce the number of lags')
else
    garg = 0.5*(df+(0:-1:1-m));
    lgg = sum(gammaln(garg));
end

% posterior.S = var.u' * var.u + prior.Sigma.scale + ...
%     prior.Phi.mean' * Ai * prior.Phi.mean + ...
%     var.B' * (var.X'*var.X) * var.B ...
%     - posterior.PhiHat' * (var.X'*var.X + Ai) * posterior.PhiHat;
% FF =  posterior.S;
% log_dnsty1 = - nobs*ny/2 * log(pi);
% log_dnsty1 = log_dnsty1 + prior.Sigma.df / 2 *log(det(prior.Sigma.scale));
% log_dnsty1 = log_dnsty1 - ny/2 * log(det(prior.Phi.cov));
% log_dnsty1 = log_dnsty1 - ny/2 * log(det( var.X'*var.X + iV )); % X'X + inv(V)
% %log_dnsty1 = log_dnsty1 - ny/2 * log(det( eye(nk) + Fv'*var.X'*var.X*Fv ));
% % % Giannone, Lenza Primiceri (2015) Appendix  A13-A14
% FF    = prior.Sigma.scale ...
%     + (posterior.PhiHat - prior.Phi.mean)'* iV * (posterior.PhiHat - prior.Phi.mean) ...
%     + var.e'* var.e;
% % FF    = eye(ny) ...
% %     + Fo'* ((posterior.PhiHat - prior.Phi.mean)'* iV * (posterior.PhiHat - prior.Phi.mean) ...
% %     + var.e'* var.e) * Fo;
% log_dnsty1 = log_dnsty1 - (nobs + prior.Sigma.df)/2 * log(det(  FF ));
% log_dnsty1 = log_dnsty1 + ggammaln(ny,(nobs + prior.Sigma.df)/2) - ggammaln(ny,prior.Sigma.df/2);
