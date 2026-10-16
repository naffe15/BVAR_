function minus_log_dnsty = bvar_opt_plr(hyperpara,y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bvar_opt_plr' computes the marginal likelihood over the hyperparameters
% of the Minnesota prior and the tightness phi of the Prior for the Long
% Run (Giannone, Lenza and Primiceri, 2019)
% Bridge function between 'bvar_ml' and the minimization routine
% 'bvar_max_hyper' (options.objective_function = 'bvar_opt_plr')

% Inputs:
% - hyperpara, [tau decay lambda mu omega phi], with one phi for all the
% rows of H or one per row; lambda and mu must be 0
% - y, data columns variables
% - lags, lag order of the VAR
% - options, as in bvar_max_hyper, plus options.priors.PLR.H (required) and
% options.priors.PLR.hyperprior (1, default: Gamma hyperprior on phi; 0 none)

% Output: minus the log marginal data density (plus the log hyperprior);
% the estimated phi are kept in [1e-5, 50]

% Gianluca Tisci, 9/26/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hyperpara(options.index_est)   = exp(hyperpara);
if isempty(options.index_fixed) == 0
    hyperpara(options.index_fixed) = exp(options.hyperpara_fidex);
end

ny = size(y,2);
if isfield(options,'priors')==0 && isfield(options,'prior')==1
    options.priors = options.prior;
end
if isfield(options,'priors')==0 || isfield(options.priors,'PLR')==0 || isfield(options.priors.PLR,'H')==0
    error('bvar_opt_plr: you must supply options.priors.PLR.H, an %dx%d invertible matrix.', ny, ny);
end
H = options.priors.PLR.H;
if size(H,1) ~= ny || size(H,2) ~= ny
    error('bvar_opt_plr: H must be %dx%d (got %dx%d).', ny, ny, size(H,1), size(H,2));
end
if rank(H) < ny
    error('bvar_opt_plr: H must be invertible (rank %d < %d).', rank(H), ny);
end
if hyperpara(3) ~= 0 || hyperpara(4) ~= 0 || any(ismember([3 4],options.index_est))
    error(['bvar_opt_plr: the PLR prior replaces the sum-of-coefficients and co-persistence ' ...
           'dummies. Set hyperpara(3:4) = 0 and leave 3 and 4 out of options.index_est.']);
end
if isfield(options,'presample')==1 && options.presample > 0
    error('bvar_opt_plr: phi tuning requires options.presample = 0.');
end

nphi = length(hyperpara) - 5;
if nphi == 1
    phi = hyperpara(6)*ones(ny,1);
elseif nphi == ny
    phi = hyperpara(6:end);
else
    error('bvar_opt_plr: hyperpara must have 6 (one phi) or %d (one phi per row of H) elements (got %d).', 5+ny, length(hyperpara));
end
if ~isreal(phi) || any(phi <= 0)
    error('bvar_opt_plr: phi must be > 0 (phi = Inf switches a row of H off).');
end
options.plr.H   = H;
options.plr.phi = phi(:);

% bounds of GLP (2019) on the estimated phi
index_phi = intersect(options.index_est, 6:length(hyperpara));
if any(hyperpara(index_phi) < 1e-5) || any(hyperpara(index_phi) > 50)
    minus_log_dnsty = 1e10;
    return
end

log_dnsty = bvar_ml(hyperpara,y,lags,options);

% Gamma hyperprior on the estimated phi, mode 1 and std 1 (GLP 2019)
plr_hyperprior = 1;
if isfield(options.priors.PLR,'hyperprior')==1
    plr_hyperprior = options.priors.PLR.hyperprior;
    if ~isequal(plr_hyperprior,0) && ~isequal(plr_hyperprior,1)
        error('bvar_opt_plr: options.priors.PLR.hyperprior must be 0 or 1.');
    end
end
if plr_hyperprior == 1
    k         = (3 + sqrt(5))/2;       % shape
    theta     = 1/sqrt(k);             % scale
    phi_est   = hyperpara(index_phi);
    log_dnsty = log_dnsty + sum((k-1)*log(phi_est) - phi_est/theta - k*log(theta) - gammaln(k));
end

minus_log_dnsty = -log_dnsty;
