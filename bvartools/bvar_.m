function [BVAR] = bvar_(y,lags,options)
%BVAR_ Bayesian VAR: posterior draws, impulse responses, forecasts, marginal likelihood.
%
%   BVAR = BVAR_(y,lags)
%   BVAR = BVAR_(y,lags,options)
%
%   Without options, a flat (Jeffreys) prior is used. All features are
%   switched on through fields of options; the main ones are listed below
%   and all are described in Chapter 3 of the Hitchhiker's Guide
%   (HitchhikerGuide_.pdf in the toolbox folder).
%
%   Inputs
%     y        T x ny (x nunits)   data, one variable per column; NaN marks
%                                  missing values (mixed frequency, ragged
%                                  edge); a third dimension stacks the units
%                                  of a panel
%     lags     scalar              number of lags (>= 1)
%     options  struct (optional), main fields:
%       Sampling and horizons
%         .K           posterior draws (default 5000)
%         .hor         IRF horizon (default 24)
%         .fhor        forecast horizon (default 12)
%       Sample and deterministic terms
%         .firstobs, .presample, .noconstant, .timetrend, .exogenous
%       Priors
%         .priors.name 'Minnesota', 'Conjugate' or 'PLR'; with
%                      .minn_prior_tau, _decay, _lambda, _mu, _omega,
%                      .max_minn_hyper and .max_plr_hyper to estimate
%                      the hyperparameters; .heterosked_weights for the
%                      Lenza-Primiceri rescaling
%       Identification (each adds its own output fields)
%         .signs, .narrative, .zeros_signs, .set, .proxy, .long_run_irf,
%         .heterosked_regimes, .hmoments, .hmoments_eig, .irf_1STD,
%         .robust_credible_regions
%       Forecasts and networks
%         .endo_path, .endo_index, .exo_index   conditional forecasts
%         .connectedness, .nethor               spillover tables
%
%   Output (main fields)
%     BVAR.Phi_draws     nk x ny x K             VAR coefficients; rows ordered
%                                               as the regressors: lags (lag-1
%                                               block first), then constant,
%                                               time trend, options.exogenous,
%                                               pandemic dummies (if present)
%     BVAR.Sigma_draws   ny x ny x K             residual covariance
%     BVAR.e_draws       residuals for each draw
%     BVAR.ir_draws      ny x hor x ny x K       recursive (Cholesky) IRFs
%     BVAR.irlr_draws, .irsign_draws, .irnarrsign_draws,
%     .irzerosign_draws, .irfset_draws, .irproxy_draws,
%     .irheterosked_draws, .irhmomsign_draws, .irhmomeig_draws
%                                                IRFs of each identification
%     BVAR.Omegas, .Omegan, .Omegaz, .Omegaset ...  accepted rotations
%     BVAR.irfset_weights                        importance weights of
%                                                options.set (use WQUANTILE)
%     BVAR.forecasts     .no_shocks, .with_shocks, .conditional
%     BVAR.logmlike      log marginal likelihood
%     BVAR.Phi_ols, .Sigma_ols, .e_ols, .InfoCrit (AIC, HQIC, BIC)
%     BVAR.Connectedness spillover tables, if requested
%
%   Output with options.exogenous_block = z (nz variables; see Guide
%   Sec. 3.1.6 and App. A.16): draws refer to the full (ny+nz) system
%     BVAR.Phi_draws     ((ny+nz)*lags+nx) x (ny+nz) x K; rows
%                        [y(-1)..y(-lags), z(-1)..z(-lags), const] (NOT
%                        ordered by lag), columns [y eqs, z eqs]
%     BVAR.Sigma_draws   (ny+nz) x (ny+nz) x K, block diagonal
%     BVAR.ir_draws      (ny+nz) x hor x (ny+nz) x K, variables and shocks
%                        ordered [y, z]
%
%   See also BVAR, CVAR_, PARSE_BVAR_OPTIONS, BVAR_MAX_HYPER, PLOT_IRFS_,
%   PLOT_FRCST_, HISTDECOMP, FEVD_.

%   Filippo Ferroni and Fabio Canova. https://github.com/naffe15/BVAR_
%
%   ROADMAP (section numbers [1]..[10] appear in the banners below)
%     [1]  options parsing        unpack opt into local variables
%     [2]  consistency checks
%     [3]  data setup and OLS     yy/XX, ydata/xdata, varols, Ridge/Lasso/EN
%     [4]  prior                  prior dummy observations and prior moments
%     [5]  posterior              posterior_ (and z-block), robust moments
%     [6]  marginal likelihood
%     [7]  posterior draws loop   Sigma, Phi -> IRFs -> forecasts ->
%                                 mixed frequency -> connectedness
%     [8]  classical and penalized estimates (OLS, Ridge, Lasso, Elastic-Net)
%     [9]  storing the Bayesian results
%     [10] nested function posterior_
%
%   PRIOR CODES (variable 'dummy', set in parse_bvar_options)
%     0 = Jeffrey / flat (default)    1 = Minnesota (dummy observations)
%     2 = conjugate MN-IW             3 = PLR (Minnesota + long-run prior)
%     4 = pandemic Minnesota
%
%   KNOWN ISSUES are tagged 'TODO(' in the code: grep for it.

if nargin < 2
    error('the bvar_ funtion needs at least two inputs: data and number of lags');
end
if lags < 1
    error('lags cannot be zero or negative');
end
%==========================================================================
% [1] OPTIONS PARSING
%   parse_bvar_options sets defaults and checks options; its output
%   struct opt is unpacked here into local variables.
%==========================================================================
if nargin > 2
    opt = parse_bvar_options(y, lags, options);
else
    opt = parse_bvar_options(y, lags);
end

% --- always set ---
ny                       = opt.ny;
nunits                   = opt.nunits;
isMatlab                 = opt.isMatlab;
K                        = opt.K;
hor                      = opt.hor;
fhor                     = opt.fhor;
nethor                   = opt.nethor;
firstobs                 = opt.firstobs;
presample                = opt.presample;
noconstant               = opt.noconstant;
timetrend                = opt.timetrend;
minn_prior_tau           = opt.minn_prior_tau;
minn_prior_decay         = opt.minn_prior_decay;
minn_prior_lambda        = opt.minn_prior_lambda;
minn_prior_mu            = opt.minn_prior_mu;
minn_prior_omega         = opt.minn_prior_omega;
long_run_irf             = opt.long_run_irf;
irf_1STD                 = opt.irf_1STD;
cfrcst_yes               = opt.cfrcst_yes;
non_explosive_           = opt.non_explosive_;
heterosked               = opt.heterosked;
signs_irf                = opt.signs_irf;
narrative_signs_irf      = opt.narrative_signs_irf;
zeros_signs_irf          = opt.zeros_signs_irf;
setrestr_irf             = opt.setrestr_irf;
proxy_irf                = opt.proxy_irf;
heterosked_irf           = opt.heterosked_irf;
hmoments_signs_irf       = opt.hmoments_signs_irf;
hmoments_eig_irf         = opt.hmoments_eig_irf;
noprint                  = opt.noprint;
nexogenous               = opt.nexogenous;
exogenous                = opt.exogenous;
cnnctdnss_               = opt.cnnctdnss_;
Ridge_                   = opt.Ridge_;
Lasso_                   = opt.Lasso_;
ElasticNet_              = opt.ElasticNet_;
set_irf                  = opt.set_irf;
robust_bayes_            = opt.robust_bayes_;
robust_credible_regions_ = opt.robust_credible_regions_;
exogenous_block          = opt.exogenous_block;
nz                       = opt.nz;
mixed_freq_on            = opt.mixed_freq_on;
dummy                    = opt.dummy;
flat                     = opt.flat;
priors                   = opt.priors;
varnames                 = opt.varnames;
y                        = opt.y;   % possibly interpolated -- see parse_bvar_options.m
pandemic_on              = opt.pandemic_on;
pandemic_h               = opt.pandemic_h;
pandemic_start           = opt.pandemic_start;
pandemic_phi             = opt.pandemic_phi;
plr                      = opt.plr;   % PLR spec ([] unless options.priors.name='PLR')
if ~isempty(plr)
    % y_bar_0 (GLP 2019, eq. 7): mean of the given initial observations.
    % Computed ONCE here so that the prior and the posterior rfvar3 calls
    % append identical PLR rows for any presample (the phi-dependent
    % normalization then cancels in posterior_int - prior_int).
    plr.ybar = mean(y(firstobs-lags : firstobs-1, :), 1);
end


% --- only present if the corresponding options.* branch fired ---
if isfield(opt,'ww'),                    ww                  = opt.ww;                    end
if isfield(opt,'K_shrinkage'),           K_shrinkage         = opt.K_shrinkage;           end
if isfield(opt,'heterosked_regimes'),    heterosked_regimes  = opt.heterosked_regimes;    end
if isfield(opt,'signs'),                 signs               = opt.signs;                 end
if isfield(opt,'narrative'),             narrative           = opt.narrative;             end
if isfield(opt,'zeros_signs'),           zeros_signs         = opt.zeros_signs;           end
if isfield(opt,'setrestr'),              setrestr            = opt.setrestr;              end
if isfield(opt,'setrestr_opt'),          setrestr_opt        = opt.setrestr_opt;          end
if isfield(opt,'f'),                     f                   = opt.f;                     end
if isfield(opt,'sr'),                    sr                  = opt.sr;                    end
if isfield(opt,'var_pos'),               var_pos             = opt.var_pos;               end
if isfield(opt,'in'),                    in                  = opt.in;                    end
if isfield(opt,'inols'),                 inols               = opt.inols;                 end
if isfield(opt,'hmoments'),              hmoments            = opt.hmoments;              end
if isfield(opt,'moment'),                moment              = opt.moment;                end
if isfield(opt,'endo_path'),             endo_path           = opt.endo_path;             end
if isfield(opt,'endo_path_index'),       endo_path_index     = opt.endo_path_index;       end
if isfield(opt,'exo_index'),             exo_index           = opt.exo_index;             end
if isfield(opt,'Omega'),                 Omega               = opt.Omega;                 end
if isfield(opt,'L'),                     L                   = opt.L;                     end
if isfield(opt,'opt_GiacomoniKitagawa'), opt_GiacomoniKitagawa = opt.opt_GiacomoniKitagawa; end
if isfield(opt,'Ridge_lambda'),          Ridge_lambda        = opt.Ridge_lambda;           end
if isfield(opt,'Lasso_lambda'),          Lasso_lambda        = opt.Lasso_lambda;           end
if isfield(opt,'ElasticNet_lambda'),     ElasticNet_lambda   = opt.ElasticNet_lambda;      end
if isfield(opt,'ElasticNet_alpha'),      ElasticNet_alpha    = opt.ElasticNet_alpha;       end
if isfield(opt,'index_nan_var'),         index_nan_var       = opt.index_nan_var;          end
if isfield(opt,'yoriginal'),             yoriginal           = opt.yoriginal;              end
if isfield(opt,'yinterpol'),             yinterpol           = opt.yinterpol;              end
if isfield(opt,'index'),                 index               = opt.index;                  end
if isfield(opt,'prior'),                 prior               = opt.prior;                  end
clear opt


%==========================================================================
% [2] CONSISTENCY CHECKS
%==========================================================================
nobs = size(y,1)-firstobs+1;
if (firstobs+ nobs-1)> size(y,1)
    fprintf('Incorrect or missing specification of the number of observations. nobs can be at most %4u\n',size(y,1)-firstobs+1);
    error('Inconsistent number of observations.')
end
if firstobs + presample + lags >= nobs  
    error('presample too large')
end
if firstobs + presample  <= lags
    error('firstobs+presample should be > # lags (for initializating the VAR)')
end
if ~isempty(plr) && mixed_freq_on == 1
    error(['PLR: missing observations are not supported. The mixed-frequency ' ...
           'Kalman path initializes from (I-B(1))\c and solves a Lyapunov ' ...
           'equation (var2ss.m, kfilternan.m), which blow up on the ' ...
           'near-unit-root draws PLR admits.']);
end
if (dummy == 1 || dummy == 3) && nexogenous > 0
    warning('I will not use exogenous variables with Minnesota');
    nexogenous = 0;
end
if mixed_freq_on == 1 && nexogenous > 0
    warning('I will not use exogenous variables with missing observations');
    nexogenous = 0;
end
if mixed_freq_on == 1 && timetrend > 0
    warning('I will not use time trend with missing observations');   
    timetrend = 0;
end
if cfrcst_yes ~= 0 && nexogenous > 0
    warning('I will not use exogenous variables with conditional forecasts');
    nexogenous = 0;
end
if size(exogenous,1) == size(y,1) && nexogenous > 0
    warning('For forecast purposes, I will assume that exo are zero out-of sample.')
    warning('To change this, include the exogenous forecasts in options.exogenous.')
    exogenous = [exogenous; zeros(fhor,nexogenous)];
end

%==========================================================================
% [3] DATA SETUP AND OLS ESTIMATES
%   Data objects:
%     idx           rows of y used in estimation (initial lags included)
%     ydata, xdata  T x ny data and T x (deterministic/exogenous) regressors,
%                   passed to rfvar3, which builds the lags itself
%     yy, XX        (T-lags) x ny and (T-lags) x nk: the regression
%                   yy = XX*Phi + E with lags built here and NO prior dummy
%                   observations; used for the residuals in the draws loop
%     varols        OLS on actual data only (no dummy observations)
%     var           OLS on actual data + prior dummy observations (posterior_)
%     nx            1 if a constant is included, 0 otherwise
%     nk            # regressors per equation (set before the draws loop)
%==========================================================================

idx = firstobs+presample-lags:firstobs+nobs-1;
nx  = 1;
if noconstant
    nx = 0;
end

% organize data as yy = XX*Phi + E (lags built here, no dummy observations)
if nunits == 1
    % EXOGENOUS BLOCK: stacked regression used ONLY to compute residuals
    % (errors = yy - XX*Phi) for each draw, not for the posterior itself.
    %   yy = [y_t , z_t]                                 (T-lags) x (ny+nz)
    %   XX = [y(-1)..y(-lags), z(-1)..z(-lags), const]   (XX2 = [z lags, const])
    % The row ordering of XX is the row ordering of the assembled Phi below.
    if exogenous_block == 1
        [yy1,XX1] = YXB_(y(idx, :),lags,[nx timetrend]);
        [yy2,XX2] = YXB_(exogenous(idx, :),lags,[nx timetrend]);
        yy = [yy1, yy2]; XX = [XX1(:,1:lags*ny), XX2];
    else 
        [yy,XX] = YXB_(y(idx, :),lags,[nx timetrend]);
    end
else % pooled units
    yy = []; XX = [];
    for nunit = 1 : nunits
        [yy_p,XX_p] = YXB_(y(idx, :, nunit),lags,[nx timetrend]);
        XX = [XX; XX_p];
        yy = [yy; yy_p];
    end
end

ydum    = [];
xdum    = [];
pbreaks = 0;
lambda  = 0;
mu      = 0;

ydata   = y(idx, :, :);
T       = size(ydata, 1);
if T-lags < lags*ny + nx
    error('Less observations than regressors: increase the # of obs or decrease the # of lags.')
end
xdata   = ones(T,nx);


if timetrend == 1
    xdata = [xdata [1-lags : T-lags]'];
end

if nexogenous > 0
    xdata = [xdata exogenous(idx,:)]; 
    XX    = [XX exogenous(idx(1)+lags : idx(end),:)];
end
% EXOGENOUS BLOCK: the lags of z enter the y-equations as additional
% 'exogenous' regressors (first lags rows are NaN = initial conditions,
% dropped by rfvar3). The same is done in posterior_ below.
if exogenous_block == 1
    xdata = [xdata,  [nan(lags,lags*nz); XX2(:,1:lags*nz)]];
end


panDum = zeros(size(y,1), pandemic_h);
for j = 1 : pandemic_h
    panDum(pandemic_start + j - 1, j) = 1;
end
xdata = [xdata, panDum(idx, :)];
XX    = [XX, panDum(idx(1)+lags : idx(end),:)]; 


% OLS estimate [NO DUMMY]:
if nunits == 1
    if heterosked == 0
        varols  = rfvar3(ydata, lags, xdata, [T; T], 0, 0);
    else
        varols  = rfvar3(ydata, lags, xdata, [T; T], 0, 0, ww);
    end
else % pooled units
    varols.y = []; varols.X = [];
    for nunit = 1 : nunits
        tmp_varols = rfvar3(ydata(:,:,nunit), lags, xdata, [T; T], 0, 0);
        varols.X = [varols.X; tmp_varols.X];
        varols.y = [varols.y; tmp_varols.y];
    end
    % Compute OLS regression and residuals of the pooled estimator
    varols = fast_ols_(varols.y,varols.X);
end
% EXOGENOUS BLOCK: OLS of the z-block VAR (z on its own lags + const).
% zvarols is NOT used for the posterior draws (those use zposterior,
% computed by posterior_(exogenous) below); zdata is used for the
% forecast initial values.
if exogenous_block == 1
    zdata   = exogenous(idx, :);
    zxdata   = ones(T,nx);
    zvarols  = rfvar3(zdata, lags, zxdata, [T; T], 0, 0);
end 

if Ridge_ == 1
    varRidge   = varols;
    varRidge.B = (eye(size(varRidge.B,1)) + Ridge_lambda*varols.xxi)\varRidge.B;
    varRidge.u = varRidge.y - varRidge.X*varRidge.B;
end
if Lasso_ == 1
    varLasso   = varols;
    for vv = 1 : ny
        varLasso.B(:,vv) = lasso(varLasso.X,varLasso.y(:,vv),'Lambda',Lasso_lambda);
    end
    varLasso.u = varLasso.y - varLasso.X*varLasso.B;
end
if ElasticNet_ == 1
    varElasticNet   = varols;
    for vv = 1 : ny
        varElasticNet.B(:,vv) = ...
            lasso(varElasticNet.X,varElasticNet.y(:,vv),'Lambda',ElasticNet_lambda,'Alpha',ElasticNet_alpha);
    end
    varElasticNet.u = varElasticNet.y - varElasticNet.X*varElasticNet.B;
end

%==========================================================================
% [4] PRIOR
%   dummy = 1, 3, 4: prior dummy observations (ydum, xdum) and the prior
%   moments obtained by OLS on presample + dummies; dummy = 0, 2: no
%   dummy observations (see Guide App. A.3-A.7).
%==========================================================================
if dummy == 0
    % JEFFREY OR UNIFORMATIVE PRIOR:
    % [priors] = jeffrey(y,lags);
    prior.name  = 'Jeffrey';
    
elseif dummy == 1
    % MINNESOTA PRIOR:
    mnprior.tight         = minn_prior_tau;
    mnprior.decay         = minn_prior_decay;
    % Use only initializations lags for the variance prior
    vprior.sig            = std(y(firstobs-lags : firstobs+presample-1,:))';
    vprior.w              = minn_prior_omega;
    lambda                = minn_prior_lambda;
    mu                    = minn_prior_mu;
    [ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior);
    % Prior density
    Tp = presample + lags;
    if nx
        xdata = xdata(1:Tp, :);
    else
        xdata = [];
    end
    varp           = rfvar3([y(firstobs-lags : firstobs+presample-1, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
    Tup            = size(varp.u, 1);
    prior.df       = Tup - ny*lags - nx - flat*(ny+1);
    prior.S        = varp.u' * varp.u;
    prior.XXi      = varp.xxi;
    prior.PhiHat   = varp.B;
    if prior.df < ny
        error('Too few degrees of freedom in the Inverse-Wishart part of prior distribution. You should increase training sample size.')
    end
    prior.minn_prior_tau    = minn_prior_tau;
    prior.minn_prior_decay  = minn_prior_decay;
    prior.minn_prior_lambda = minn_prior_lambda;
    prior.minn_prior_mu     = minn_prior_mu;
    prior.minn_prior_omega  = minn_prior_omega;
       
elseif dummy == 2
    % MN-IW
    prior.name  = 'MultivariateNormal-InverseWishart';    

elseif dummy == 3
    % PLR PRIOR: Minnesota short-run dummies + PLR long-run dummies
    mnprior.tight         = minn_prior_tau;
    mnprior.decay         = minn_prior_decay;
    % Use only initializations lags for the variance prior
    vprior.sig            = std(y(firstobs-lags : firstobs+presample-1,:))';
    vprior.w              = minn_prior_omega;
    lambda                = minn_prior_lambda;
    mu                    = minn_prior_mu;
    [ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior);
    % Prior density
    Tp = presample + lags;
    if nx
        xdata = xdata(1:Tp, :);
    else
        xdata = [];
    end
    varp           = rfvar3([y(firstobs-lags : firstobs+presample-1, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu, [], plr);
    Tup            = size(varp.u, 1);
    prior.df       = Tup - ny*lags - nx - flat*(ny+1);
    prior.S        = varp.u' * varp.u;
    prior.XXi      = varp.xxi;
    prior.PhiHat   = varp.B;
    if prior.df < ny
        error('PLR: too few degrees of freedom in the Inverse-Wishart part of the prior. Tighten the Minnesota short-run or increase the training sample size.')
    end
    prior.minn_prior_tau    = minn_prior_tau;
    prior.minn_prior_decay  = minn_prior_decay;
    prior.minn_prior_lambda = minn_prior_lambda;
    prior.minn_prior_mu     = minn_prior_mu;
    prior.minn_prior_omega  = minn_prior_omega;
    prior.plr               = plr;

elseif dummy == 4
    mnprior.tight         = minn_prior_tau;
    mnprior.decay         = minn_prior_decay;
    vprior.sig            = std(y(firstobs-lags : firstobs+presample-1,:))';
    vprior.w              = minn_prior_omega;
    lambda                = minn_prior_lambda;
    mu                    = minn_prior_mu;
    pandemic.h            = pandemic_h;
    pandemic.phi          = pandemic_phi;
    [ydum, xdum, pbreaks] = pandemic_varprior(ny, nx, lags, mnprior, vprior, pandemic);
    % Prior density
    Tp = presample + lags;
    if nx
        xdata = xdata(1:Tp, :);
    else
        xdata = [];
    end
    varp           = rfvar3([y(firstobs-lags : firstobs+presample-1, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
    Tup            = size(varp.u, 1);
    prior.df       = Tup - ny*lags - nx - flat*(ny+1) - pandemic_h;
    prior.S        = varp.u' * varp.u;
    prior.XXi      = varp.xxi;
    prior.PhiHat   = varp.B;
    if prior.df < ny
        error('Too few degrees of freedom in the Inverse-Wishart part of prior distribution. You should increase training sample size.')
    end
    prior.minn_prior_tau    = minn_prior_tau;
    prior.minn_prior_decay  = minn_prior_decay;
    prior.minn_prior_lambda = minn_prior_lambda;
    prior.minn_prior_mu     = minn_prior_mu;
    prior.minn_prior_omega  = minn_prior_omega;
    prior.pandemic_h        = pandemic_h;
    prior.pandemic_phi      = pandemic_phi;


end

%==========================================================================
% [5] POSTERIOR
%   posterior_ re-runs OLS on actual data + prior dummy observations and
%   returns the Normal-Inverse-Wishart posterior moments (see [10]).
%==========================================================================
[posterior,var] = posterior_(y);

% EXOGENOUS BLOCK (Guide Sec. 3.1.6 and App. A.16; ny = n, nz = q there):
% the block-exogenous VAR is treated as TWO INDEPENDENT
% VARs, each with its own Normal-Inverse-Wishart posterior:
%   (a) y-block:  y_t on [y lags, const, z lags]  -> posterior (computed above)
%   (b) z-block:  z_t on [z lags, const]          -> zposterior (computed here)
% The factorisation holds because (i) y lags do not enter the z equations
% and (ii) Cov(u_t, e_t) = 0, so the likelihood (and, with the flat prior,
% the posterior) splits into a y part and a z part.
% To reuse posterior_ for (b) we temporarily switch the flag off (so that
% posterior_ does not append the z lags to xdata) and set ny = nz; both are
% restored right after.
% TODO(exo_block): the Minnesota/conjugate dummies (ydum, xdum, prior.Phi)
% are built for the y-block only: with the exogenous block only the flat
% prior is consistent at the moment.
% TODO(exo_block): with the flat prior the y-block posterior.df is
% Tu-ny*lags-nx+ny+1, i.e. it ignores the nz*lags z-lag regressors
% (should be T-k_y+n+1 with k_y = (n+q)*lags+1, Guide App. A.16).
% TODO(exo_block): robust_bayes_ > 0 does not define ySigma/zSigma.
if exogenous_block == 1
    exogenous_block = 0;
    ny0 = ny;
    ydata0 = ydata;
    ny = nz;
    [zposterior, zvar] = posterior_(exogenous);
    exogenous_block = 1;
    ny = ny0;
    ydata = ydata0;
end

% compute the second and fourth moments robust to error distribution misspecification
% (Petrova 2022; Guide App. A.8)  
if robust_bayes_ > 0 
    
    % Kurtosis
    Dpl       = duplicationmatrix(ny);
    Dplus     = pinv(Dpl);
    I_ny      = vech(eye(ny));
    Sig_      = 1/nobs * varols.u' * varols.u;
    vech_Sig_ = vech(Sig_);

    % whiten the reduced form errors
    SigChol   = chol(Sig_,'lower');
    iSigChol  = inv(SigChol);
    white_u   = varols.u * iSigChol';     
    % compute fourth centered moments
    Khat_     = fourthmom(white_u);
    % shrink towards the normal multivariate fourth mom
    % Petrova (2022), equation (17) at page 161 (TODO(ref): check)
    K1_       = nobs/(nobs+K_shrinkage) * (Khat_ - I_ny * I_ny') ; 
    K2_       = K_shrinkage/(nobs+K_shrinkage) * ...
        Dplus *(eye(ny^2) + commutationmatrix(ny)) * Dplus'; 
    Kstar_     = K1_ + K2_;
    % Petrova (2022), equation (8) of remark 3 at page 158 (TODO(ref): check)
    Left  = Dplus * kron(SigChol,SigChol) * Dpl;
    Rigth = Left';
    
    vech_Sig_cov_            = 1/nobs * Left * (Kstar_) * Rigth;   
    vech_Sig_cov_lower_chol  = chol(vech_Sig_cov_,'lower');
    
    if robust_bayes_ > 1 
        % Skewness
        ivech_Sig_cov_ = pinv(vech_Sig_cov_);
        Sstar_         = nobs /(nobs+K_shrinkage) * thirdmom(varols.u);
        mu_cov_        = 1/nobs * Sstar_ * ivech_Sig_cov_ * Sstar_';
        mu_cov_chol    = chol(Sig_ - mu_cov_,'lower');
    end
end

%==========================================================================
% [6] LOG MARGINAL LIKELIHOOD (log data density)
%   dummy = 1, 3, 4: log ratio of the matrix-t integrating constants of
%     posterior and prior (matrictint), minus the Gaussian constant;
%   dummy = 0: posterior constant only; dummy = 2: mniw_log_dnsty.
%   TODO(exo_block): only the y-block enters log_dnsty.
%==========================================================================

try
    posterior_int = matrictint(posterior.S, posterior.df, posterior.XXi);
catch
    warning('I could not compute the marginal likelihood');
    posterior_int = nan;
end
        
if dummy == 0 % jeffrey prior
    lik_nobs  = posterior.df;
    log_dnsty = posterior_int - 0.5*ny*lik_nobs*log(2*pi);

elseif dummy == 1 % only for minnesota dummy
    prior_int = matrictint(prior.S, prior.df, prior.XXi);
    lik_nobs  = posterior.df - prior.df;
    log_dnsty = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);

elseif dummy == 2 % conjugate MN-IW prior
    log_dnsty = mniw_log_dnsty(prior,posterior,varols);

elseif dummy == 3 % PLR: same conjugate NIW marginal as the Minnesota dummy prior
    prior_int = matrictint(prior.S, prior.df, prior.XXi);
    lik_nobs  = posterior.df - prior.df;
    log_dnsty = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);

elseif dummy == 4
    prior_int = matrictint(prior.S, prior.df, prior.XXi);
    lik_nobs  = posterior.df - prior.df;
    log_dnsty = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);  
    
end


%==========================================================================
% [7] POSTERIOR DRAWS
%   Preallocation, then for each draw d = 1..K: [7a] draw (Sigma, Phi),
%   [7b] IRFs, [7c] forecasts, [7d] mixed frequency, [7e] connectedness.
%   The last dimension of every *_draws array indexes the draw.
%==========================================================================

% Preallocation of memory
% Matrices for collecting draws from Posterior Density
% Last dimension corresponds to a specific draw

% lags, const, trend, exogenous, pandemic dummies
Phi_draws = zeros(ny*lags+nx+timetrend+nexogenous+pandemic_h, ny, K); 
Sigma_draws   = zeros(ny,ny,K);                     % Shocks Covariance
Sigma_lower_chol_draw = zeros(ny,ny,K);             % Cholesky factor of Sigma 
ir_draws      = zeros(ny,hor,ny,K);                 % variable, horizon, shock and draws - Cholesky IRF
irlr_draws    = zeros(ny,hor,ny,K);                 % variable, horizon, shock and draws - Long Run IRF
Qlr_draws     = zeros(ny,ny,K);                     % long run impact matrix
e_draws       = zeros(size(yy,1), ny,K);                  % residuals
% in-sample forecast errors are stored only if the array has < 1e8
% elements (about 0.8 GB in double precision)
if size(yy,1)*ny*fhor*K < 1e8
    flag_fe  = 1;
    fe_draws = zeros(size(yy,1), ny,fhor,K);  % forecast errors 
else 
    flag_fe  = 0; 
end
yhatfut_no_shocks         = NaN(fhor, ny, K, nunits);   % forecasts without shocks
yhatfut_with_shocks       = NaN(fhor, ny, K, nunits);   % forecasts with simulated shocks
yhatfut_cfrcst            = NaN(fhor, ny, K, nunits);   % forecast conditional on endogenous path
logL                      = NaN(K,1);
% EXOGENOUS BLOCK: all draws are stored for the full (ny+nz) system.
if exogenous_block == 1
    Phi_draws     = zeros((ny+nz)*lags+nx, (ny+nz), K);
    Sigma_draws   = zeros((ny+nz),(ny+nz),K);
    Sigma_lower_chol_draw = zeros((ny+nz),(ny+nz),K);             % Cholesky factor of Sigma 
    ir_draws      = zeros((ny+nz),hor,(ny+nz),K);
    irlr_draws    = zeros((ny+nz),hor,(ny+nz),K);
    Qlr_draws     = zeros((ny+nz),(ny+nz),K);
    e_draws       = zeros(size(yy,1),(ny+nz),K);
    fe_draws      = zeros(size(yy,1),(ny+nz),fhor,K);
    yhatfut_no_shocks         = NaN(fhor, (ny+nz), K, nunits);
    yhatfut_with_shocks       = NaN(fhor, (ny+nz), K, nunits);
    yhatfut_cfrcst            = NaN(fhor, (ny+nz), K, nunits);
end 
if signs_irf == 1
    irsign_draws = ir_draws;
    Omega_draws  = Sigma_draws;
    OmegaEmpty   = zeros(K, 1);
    if robust_credible_regions_
        irsign_lower_draws = ir_draws;
        irsign_upper_draws = ir_draws;
    end
end
if nexogenous > 0 
    irx_draws = zeros(ny,hor,nexogenous,K);  
end
if narrative_signs_irf == 1
    irnarrsign_draws = ir_draws;
    Omegan_draws     = Sigma_draws;
    OmegaEmpty       = zeros(K, 1);
    if robust_credible_regions_
        irnarrsign_lower_draws = ir_draws;
        irnarrsign_upper_draws = ir_draws;
    end
end
if zeros_signs_irf == 1
    irzerosign_draws   = ir_draws;
    Omegaz_draws       = Sigma_draws;
    OmegaEmpty         = zeros(K, 1);
    if robust_credible_regions_
        irzerosign_lower_draws = ir_draws;
        irzerosign_upper_draws = ir_draws;
    end
end
if setrestr_irf == 1
    irfset_draws       = ir_draws;
    Omegaset_draws     = Sigma_draws;
    irfset_weights     = zeros(K, 1);
    OmegaEmpty         = zeros(K, 1);
    if robust_credible_regions_
        irfset_lower_draws = ir_draws;
        irfset_upper_draws = ir_draws;
    end
end
if proxy_irf == 1
    irproxy_draws = ir_draws;
    Omegap_draws  = Sigma_draws;
end
if heterosked_irf == 1
   irheterosked_draws = ir_draws;
   Omegah_draws       = Sigma_draws;
end
if hmoments_signs_irf == 1
    irhmomsign_draws = ir_draws;
    Omegam_draws     = Sigma_draws;
    OmegaEmpty         = zeros(K, 1);
    if robust_credible_regions_
        irhmomsign_lower_draws = ir_draws;
        irhmomsign_upper_draws = ir_draws;
    end
end
if hmoments_eig_irf == 1
    irhmomeig_draws = ir_draws;
    Omegae_draws     = Sigma_draws;
end
if mixed_freq_on
    yfill = nan(size(y,1),ny,K);
    yfilt = nan(size(y,1),ny,K);    
    logL  = nan(K,1);
end
if cnnctdnss_ == 1
    CnndtnssIndex         = nan(K,1);
    CnndtnssFromAlltoUnit = nan(ny,K);
    CnndtnssFromUnitToAll = nan(ny,K);
    Ctheta                = nan(ny,ny,K);
end

% Settings for the forecasts
forecast_data.xdata       = ones(fhor, nx);
if timetrend
    forecast_data.xdata = [forecast_data.xdata (T-lags+1 : T-lags+fhor)'];
end
if nexogenous>0
    forecast_data.xdata = [forecast_data.xdata exogenous(T-lags+1 : T-lags+fhor,:)];
end
forecast_data.xdata = [forecast_data.xdata, zeros(fhor, pandemic_h)];

% EXOGENOUS BLOCK: forecasts are run on the full system, so the initial
% values stack the last lags observations of [y, z].
if exogenous_block == 1
    forecast_data.initval     = [ydata(end-lags+1:end, :), zdata(end-lags+1:end, :)];
else
    forecast_data.initval     = ydata(end-lags+1:end, :, :);
end

% Settings for the MFVAR
if mixed_freq_on == 1
    KFoptions.index   = index;
    KFoptions.noprint = noprint;
end


% Cholesky of inv(posterior.S), input to rand_inverse_wishart; fallbacks if
% S is numerically not positive definite: LDL, then a 1e-05 diagonal ridge.
try
    S_inv_upper_chol    = chol(inv(posterior.S));
catch
    warning('POSTERIOR MEAN of SIGMA IS ILL-BEHAVED (NON-POSITIVE DEFINITE)')
    try
        warning('I will try with the LDL decomposition')
        iS               = inv(posterior.S);
        [Left, D, Right]        = ldl(iS);
        S_inv_upper_chol = sqrt(D) * Left' * Right';
    catch
        warning('I will add +1e-05 to the diagonal')
        S_inv_upper_chol    = chol(inv(posterior.S + 1e-05*eye(ny)));
    end
end

XXi_lower_chol      = chol(posterior.XXi)';
% EXOGENOUS BLOCK: Cholesky factors of the z-block posterior, computed once.
% ZZi_lower_chol: lower chol of (X_z'X_z)^{-1} of the z-block regression
% zS_inv_upper_chol: upper chol of inv(S_z), input to rand_inverse_wishart
if exogenous_block == 1
    ZZi_lower_chol      = chol(zposterior.XXi)';
    zS_inv_upper_chol    = chol(inv(zposterior.S));
end

% EXOGENOUS BLOCK: nk = # regressors in each y-equation
% (ny*lags y lags + nx const + timetrend + nz*lags z lags).
if exogenous_block == 1
    nk = (ny+nz)*lags+nx+timetrend;
else
    nk = ny*lags+nx+timetrend + nexogenous;
end


nk = nk + pandemic_h;


% Declaration of the companion matrix
Companion_matrix = diag(ones(ny*(lags-1),1),-ny);
% p         counts draws rejected by the stability check (not reported)
% dd        = 1 once the current draw is accepted; reset at the end of the loop
% nan_count unused (left over from debugging)
p  = 0;
dd = 0;
nan_count = 0;
waitbar_yes = 0;
% progress bar only for more than 99 draws
if K > 99
    waitbar_yes = 1;
    wb = waitbar(0, 'Generating draws from the Posterior Distribution');
end

for  d =  1 : K
    
    while dd == 0 
        %======================================================================
        % [7a] Draw (Sigma, Phi) from the posterior (Guide App. A.3-A.8)
        % Step 1: draw Sigma
        %   robust_bayes_ > 0: vech(Sigma) from its asymptotic normal
        %     distribution (Petrova 2022, Guide App. A.8), redrawn until PD
        %   otherwise: Sigma ~ IW(posterior.S, posterior.df)
        if robust_bayes_ > 0 % robust to kurtosis 
            tec = 0;
            while tec == 0
                Sig0  = randn((ny+1)*ny/2, 1);
                Sig1  = vech_Sig_cov_lower_chol * Sig0;
                Sig2  = vech_Sig_ + Sig1;
                Sigma = ivech(Sig2);
                [Sigma_lower_chol,flag] = chol(Sigma,'lower');
                if flag == 0
                    tec = 1;
                end
            end
        else  
            % EXOGENOUS BLOCK [App. A.16, drawing the disassembled blocks]:
            % draw the two covariance matrices independently
            %    ySigma ~ IW(S_y, df_y) (Sigma in the Guide), zSigma ~ IW(S_z, df_z) (Sigma_e)
            % and stack them block-diagonally (Cov(u_t,e_t) = 0 by assumption).
            % The block Cholesky factors are needed for the coefficient draws below.
            if exogenous_block == 1
                 ySigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);  
                 zSigma = rand_inverse_wishart(nz, zposterior.df, zS_inv_upper_chol);
                 Sigma = [ySigma, zeros(ny,nz); zeros(nz,ny), zSigma];
                 Sigma_lower_chol = chol(Sigma)';
                 ySigma_lower_chol = chol(ySigma)';
                 zSigma_lower_chol = chol(zSigma)';
            else
                Sigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);
                Sigma_lower_chol = chol(Sigma)';
            end
        end 
        % Step 2: draw from the AR coefficients
        % Step 2: draw Phi | Sigma ~ MN(posterior.PhiHat, posterior.XXi, Sigma)
        if exogenous_block == 1
            % Draw Phi from its matrix-normal conditional posterior without
            % ever forming the full kron(Sigma,XXi) matrix: for Z iid N(0,1)
            % reshaped to nk x ny, XXi_lower_chol*Z*Sigma_lower_chol' has the
            % same distribution as reshape(kron(Sigma_lower_chol,XXi_lower_chol)*vec(Z), nk, ny),
            % via vec(A*Z*B') = kron(B,A)*vec(Z). Same draw, far cheaper for
            % large ny*nk since the (ny*nk)x(ny*nk) Kronecker matrix is never built.
            % EXOGENOUS BLOCK [App. A.16, drawing the disassembled blocks]:
            % conditional on its own Sigma, draw each block's
            % coefficients from its matrix-normal posterior
            %   (a) y-block: yPhi ~ MN(PhiHat_y, (X_y'X_y)^{-1}, ySigma),  nk x ny
            %       rows of yPhi: [y(-1)..y(-lags); const (; trend); z(-1)..z(-lags)]
            Phi1 = randn(nk * ny, 1);
            Z    = reshape(Phi1, nk, ny);
            Phi3 = XXi_lower_chol * Z * ySigma_lower_chol';
            yPhi  = Phi3 + posterior.PhiHat;
            %   (b) z-block: zPhi ~ MN(PhiHat_z, (X_z'X_z)^{-1}, zSigma), (nz*lags+nx+trend) x nz
            %       rows of zPhi: [z(-1)..z(-lags); const (; trend)]
            zPhi1 = randn( (nz*lags+nx+timetrend) * nz, 1);
            Zz    = reshape(zPhi1, (nz*lags+nx+timetrend), nz);
            zPhi3 = ZZi_lower_chol * Zz * zSigma_lower_chol';
            zPhi  = zPhi3 + zposterior.PhiHat;

            % EXOGENOUS BLOCK [App. A.16, assembling the stacked form]:
            % disassemble the two draws and reassemble them
            % into the stacked coefficient matrix of the full system
            %                  y eqs (ny)   z eqs (nz)
            %   y lags     [   yPhiy        zPhiy = 0   ]   <- block exogeneity restriction
            %   z lags     [   yPhiz        zPhiz       ]
            %   const/tr   [   yPhic        zPhic       ]
            % Row ordering = columns of XX (all y lags first, then all z lags), i.e.
            % NOT lag-interleaved: this is what is stored in Phi_draws and used for
            % residuals and forecasts. The IRF section reorders it lag by lag.
            % TODO(exo_block): the slicing of yPhi assumes no options.exogenous and no pandemic
            % dummies (they would sit between const/trend and the z lags).
            yPhiy = yPhi(1:ny*lags,:); % coefficients of lag endogenous on endogenous
            yPhiz = yPhi(ny*lags+nx+timetrend+1:end,:); % coefficient of lag exogenous on endogenous
            yPhic = yPhi(ny*lags+1 : ny*lags+nx+timetrend, :); % constant and time trend of endogenous
            zPhiz = zPhi(1:nz*lags,:); % coefficient of lag exogenous on exogenous
            zPhiy = zeros(ny*lags, nz); % coefficient of lag endogenous on exogenous
            zPhic = zPhi(nz*lags+1 : nz*lags+nx+timetrend, :);% constant and time trend of exogenous
            % 
            Phi = [yPhiy, zPhiy;
                    yPhiz, zPhiz;
                    yPhic, zPhic];
        else
            % See comment above: avoids forming kron(Sigma_lower_chol,XXi_lower_chol).
            Phi1 = randn(nk * ny, 1);
            Z    = reshape(Phi1, nk, ny);
            Phi3 = XXi_lower_chol * Z * Sigma_lower_chol';
            Phi  = Phi3 + posterior.PhiHat;
        end
        
        % With skewness, the intercept is redrawn from its asymptotic
        % distribution conditional on Sigma (Petrova 2022, Guide App. A.8).
        if robust_bayes_ > 1 % robust to kurtosis and skewness
            mu_rob           = posterior.PhiHat(ny*lags+1,:) + (Sstar_* ivech_Sig_cov_ * (Sig2  - vech_Sig_))';
            % as if assuming no lags (only costant) (X'X) = T
            % mu_cov_chol      = chol(Sigma - mu_cov_,'lower');            
            Phi(ny*lags+1,:) = mu_rob + (mu_cov_chol*randn(ny, 1))';
        end
        
        % Step 3: stability check (only if non_explosive_ = 1). The companion
        % eigenvalues in 'test' are reused in [7d] (mixed frequency), so the
        % eig is computed only when one of the two needs it: it costs
        % O((ny*lags)^3) per draw.
        if exogenous_block == 0
            if non_explosive_ == 1 || mixed_freq_on
                Companion_matrix(1:ny,:) = Phi(1:(ny)*lags,:)';
                test = (abs(eig(Companion_matrix)));
            end
            if non_explosive_ == 1
                % Checking the eigenvalues of the companion matrix (on or inside the
                % unit circle)
                % tolerance: roots up to 1.01 in modulus are accepted
                if any(test>1.01)
                    p = p+1;
                else
                    dd = 1;
                end
            else
                dd = 1;
            end
        else
            % TODO(exo_block): no stability check, every draw is accepted.
            dd = 1;
        end
    end
    
    % store the draws
    Phi_draws(:,:,d)   = Phi;
    Sigma_draws(:,:,d) = Sigma;
    errors             = yy - XX * Phi;
    e_draws(:,:,d)     = errors;
    Sigma_lower_chol_draw(:,:,d) = Sigma_lower_chol;
    % EXOGENOUS BLOCK: from here to the end of the loop the code works on the
    % full system, so ny is temporarily set to ny+nz (restored at the end of
    % the iteration). TODO(exo_block): fe_ is fed Phi(1:ny*lags,1:ny) in the stacked (non-
    % interleaved) ordering, while var2ma inside fe_ expects lag-1 block first:
    % fe_draws are not correct for lags > 1 (to be fixed).
    if exogenous_block == 1
        ny0 = ny;
        ny = ny + nz;
    end 
    if flag_fe
        fe_draws(:,:,:,d)  = fe_(fhor,errors,Phi(1 : ny*lags, 1 : ny));
    end
    
    %======================================================================
    % [7b] IRFs
    % EXOGENOUS BLOCK [App. A.16, from the stacked to the companion form]:
    % iresponse needs the companion ordering
    % [y(-1), z(-1), y(-2), z(-2), ...]. Phi_tmp keeps the stacked ordering
    % (restored before the forecasts); Phi0 is the lag-interleaved version.
    % With block-diagonal Sigma and zPhiy = 0, z does not respond to y shocks.
    if exogenous_block
        % reorder regressors from y(-1), y(-2), ..., z(-1), z(-2), ...
        % to  y(-1), z(-1), y(-2), z(-2) ...
        Phi_tmp = Phi; Phi0 = zeros(size(Phi));
        index = 0;
        for ell = 1 : lags
            Phi0(index+ 1 : index + ny0, :) = Phi_tmp(ny0*(ell-1)+1 : ny0*ell,:);
            Phi0(index + ny0 + 1 : index + ny0 + nz, :) = ... 
            Phi_tmp(ny0*lags + (ell-1)*nz + 1 : ny0*lags + (ell)*nz, :); 
            
            index = index + ny; 
        end
        Phi0(ny*lags+1:end, :) = Phi(ny*lags+1:end, :);
        Phi = Phi0;
    end
    % Compute the impulse response functions
    % with cholesky
    if irf_1STD == 1
        % one STD increase
        ir_draws(:,:,:,d)      = iresponse(Phi(1 : ny*lags, 1 : ny),Sigma,hor,eye(ny));
    else
        % one percent increase
        ir_draws(:,:,:,d)      = iresponse(Phi(1 : ny*lags, 1 : ny),Sigma,hor,eye(ny),0);
    end
    if nexogenous > 0
        Phi1         = Phi(1 : ny*lags + nx + timetrend, 1 : ny);
        Exo1         = zeros(ny,nexogenous);
        Exo1(:,1:nexogenous) = Phi(ny*lags + nx + timetrend +1 ....
                            : ny*lags + nx + timetrend + nexogenous, 1 : ny)';
        Sig1         = eye(ny);
        % unitary increase        
        irx_draws(:,:,:,d) = iresponse(Phi1,Sig1,hor,Exo1);        
    end
    % define the identity rotation
    % NB: each identification scheme below overwrites Omega; the Omega of the
    % LAST active scheme is the one used in [7c] by cforecasts2
    % (cfrcst_yes = 2) and in [7e] by connectedness (cnnctdnss_ = 2).
    Omega = eye(ny);
    
    % with long run restrictions
    if long_run_irf == 1
        [irlr,Qlr]             = iresponse_longrun(Phi(1 : ny*lags, 1 : ny),Sigma,hor,lags);
        irlr_draws(:,:,:,d)    = irlr;
        Qlr_draws(:,:,d)       = Qlr;
        Omega                  = Qlr;
    end
    % with sign restrictions
    if signs_irf == 1
        [irsign,Omega]         = iresponse_sign(Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs);
        irsign_draws(:,:,:,d)  = irsign;
        Omega_draws(:,:,d)     = Omega;
        OmegaEmpty(d)          = any(any(isnan(Omega)));
        % Approximating Giacomini-Kitagawa robust credible sets
        if robust_credible_regions_            
            if OmegaEmpty(d) == 0          
                irsign0 = nan(ny,hor,ny,L);
                for ell = 1 : L
                    [irsign0(:,:,:,ell), ~] = iresponse_sign(Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs);
                end
                irsign_upper_draws(:,:,:,d) = max(irsign0,[],4,"omitnan");
                irsign_lower_draws(:,:,:,d) = min(irsign0,[],4,"omitnan");
            end
        end
    end
    % with narrative and sign restrictions
    if narrative_signs_irf == 1
        [irnarrsign,Omega]         = iresponse_sign_narrative(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,narrative);
        irnarrsign_draws(:,:,:,d)  = irnarrsign;
        Omegan_draws(:,:,d)        = Omega;
        OmegaEmpty(d)          = any(any(isnan(Omega)));
        % Approximating Giacomini-Kitagawa robust credible sets
        if robust_credible_regions_
            if OmegaEmpty(d) == 0
                irsign0 = nan(ny,hor,ny,L);
                for ell = 1 : L
                    [irsign0(:,:,:,ell), ~] = iresponse_sign_narrative(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,narrative);
                end
                irnarrsign_upper_draws(:,:,:,d) = max(irsign0,[],4,"omitnan");
                irnarrsign_lower_draws(:,:,:,d) = min(irsign0,[],4,"omitnan");
            end
        end

    end
    % with zeros and sign restrictions
    if zeros_signs_irf == 1
        [irzerosign,Omega]          = iresponse_zeros_signs(Phi,Sigma,hor,lags,var_pos,f,sr);
        irzerosign_draws(:,:,:,d)   = irzerosign;
        Omegaz_draws(:,:,d)         = Omega;
        OmegaEmpty(d)          = any(any(isnan(Omega)));
        % Approximating Giacomini-Kitagawa robust credible sets
        if robust_credible_regions_
            if OmegaEmpty(d) == 0
                irsign0 = nan(ny,hor,ny,L);
                for ell = 1 : L
                    [irsign0(:,:,:,ell), ~] = iresponse_zeros_signs(Phi,Sigma,hor,lags,var_pos,f,sr);
                end
                irzerosign_upper_draws(:,:,:,d) = max(irsign0,[],4,"omitnan");
                irzerosign_lower_draws(:,:,:,d) = min(irsign0,[],4,"omitnan");
            end
        end

    end
    % with options.set: sign, zero, narrative and dominance restrictions
    if setrestr_irf == 1
        [irfset,Omega,irfset_weights(d)] = iresponse_set(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,setrestr,[],setrestr_opt);
        irfset_draws(:,:,:,d) = irfset;
        Omegaset_draws(:,:,d) = Omega;
        OmegaEmpty(d)         = any(any(isnan(Omega)));
        % Approximating Giacomini-Kitagawa robust credible sets
        if robust_credible_regions_ && OmegaEmpty(d) == 0
            optL       = setrestr_opt;
            optL.draws = L;
            irsign0    = iresponse_set(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,setrestr,[],optL);
            irfset_upper_draws(:,:,:,d) = max(irsign0,[],4);   % NaN ignored
            irfset_lower_draws(:,:,:,d) = min(irsign0,[],4);
        end
    end
    % with proxy
    if proxy_irf == 1
        in.res                  = e_draws(:,:,d);
        in.Phi                  = Phi_draws(:,:,d)  ;
        in.Sigma                = Sigma;
        tmp_                    = iresponse_proxy(in);
        irproxy_draws(:,:,1,d)  = tmp_.irs';        
        omega11 = tmp_.irs(1,1);
        omega21 = tmp_.irs(1,2:end);
        omega22 = chol(Sigma(2:end,2:end)-omega21'*omega21)' * generateQ(ny-1);
        iomega22 = inv(omega22');
        omega12o  = (Sigma(1,2:end)- omega11*omega21 )*iomega22;
        norm     = sqrt(      omega12o*omega12o');
        omega12 = omega12o/norm*sqrt(Sigma(1,1)-omega11^2); 
        icholSig = inv(Sigma_lower_chol);
        Omega = icholSig*[[omega11; omega21'] , [omega12; omega22]];
        Omegap_draws(:,:,d) = Omega;
        clear tmp_
    end
    % with heteroskedasticity 
    if heterosked_irf == 1
        [irheterosked,Omega]       = iresponse_heterosked(Phi(1 : ny*lags, 1 : ny),errors,hor,heterosked_regimes);
        irheterosked_draws(:,:,:,d) = irheterosked;
        Omegah_draws(:,:,d)         = Omega;
    end
    % with higher-moments and sign restrictions
    if hmoments_signs_irf == 1
        [irhmomsign,Omega]         = iresponse_sign_hmoments(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,hmoments,f);
        irhmomsign_draws(:,:,:,d)  = irhmomsign;
        Omegam_draws(:,:,d)        = Omega;
        OmegaEmpty(d)          = any(any(isnan(Omega)));
        if robust_credible_regions_            
            if OmegaEmpty(d) == 0          
                irsign0 = nan(ny,hor,ny,L);
                parfor ell = 1 : L
                    [irsign0(:,:,:,ell), ~] = iresponse_sign_hmoments(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,hmoments,f);
                end
                irhmomsign_upper_draws(:,:,:,d) = max(irsign0,[],4,"omitnan");
                irhmomsign_lower_draws(:,:,:,d) = min(irsign0,[],4,"omitnan");
            end
        end
    end
    % with higher-moments eigenvalue decomposition
    if hmoments_eig_irf == 1
        [irhmomeig,Omega]         = iresponse_sign_hmoments_eig(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,moment);
        irhmomeig_draws(:,:,:,d)  = irhmomeig;
        Omegae_draws(:,:,d)        = Omega;
    end
    
    %======================================================================
    % [7c] Forecasts (unconditional and conditional)
    % EXOGENOUS BLOCK: forecasts use the stacked (non-interleaved) Phi,
    % consistent with forecast_data.initval = [y, z].
    if exogenous_block
        Phi = Phi_tmp;
    end
    % compute the out of sample forecast (unconditional)
    [frcst_no_shock,frcsts_with_shocks] = forecasts(forecast_data,Phi,Sigma,fhor,lags);
    yhatfut_no_shocks(:,:,d,:)            = frcst_no_shock;
    yhatfut_with_shocks(:,:,d,:)          = frcsts_with_shocks;

    if cfrcst_yes == 1
        % Forecast conditional on the path of an endo var using all shocks
        [sims_with_endopath,EPS(:,:,d)] = ...
            cforecasts(endo_path,endo_path_index,forecast_data,Phi,Sigma);
        yhatfut_cfrcst(:,:,d) = sims_with_endopath;
        
    elseif cfrcst_yes == 2
        % Forecast conditional on the path of an endo var using only a
        % subset of shocks. notice that the # of endo and # exo must coincide
        % Omega is the structural orthonormal matrix
        [sims_with_endopath,EPS(:,:,d)] = ...
            cforecasts2(endo_path,endo_path_index,exo_index,forecast_data,Phi,Sigma,Omega);
        yhatfut_cfrcst(:,:,d) = sims_with_endopath;
    end
    
    
    %======================================================================
    % [7d] Mixed frequency: Kalman filter/smoother, then posterior update
    if mixed_freq_on
        % If the draw has a root outside the unit circle, the Kalman filter
        % is initialized differently (KFoptions.initialCond = 1).
        % TODO(mixed_freq): 'test' is NOT recomputed here; it is the one from
        % the stability check of [7a] (same draw), which is skipped with the
        % exogenous block.
        % initialCond is set draw by draw (it was never reset to 0 before).
        if any(test>1.0000000000001)
            KFoptions.initialCond = 1;
        else
            KFoptions.initialCond = 0;
        end
        % Forward Kalman Filter and Smoother
        [logL,KFout]     = kfilternan(Phi,Sigma,yoriginal,KFoptions);
        yfill(:,:,d)     = KFout.smoothSt_plus_ss(:,KFout.index_var);
        yfilt(:,:,d)     = KFout.filteredSt_plus_ss(:,KFout.index_var);
        % recompute the posterior with smoothed data
        % (prior is fixed -- presample assumed complete, not Kalman-filled)
        [posterior1]     = posterior_(yfill(:,:,d));
        S_inv_upper_chol = chol(inv(posterior1.S));
        XXi_lower_chol   = chol(posterior1.XXi)';
        posterior.PhiHat = posterior1.PhiHat;
        forecast_data.initval = yfill(end-lags+1:end, :, d);
        logL(d) = KFout.logL;
    end
    
    %======================================================================
    % [7e] Connectedness
    if cnnctdnss_ > 0 
        if cnnctdnss_ == 2 
            [C] = connectedness(Phi(1 : ny*lags, 1 : ny),Sigma,nethor, Sigma_lower_chol * Omega );
        elseif cnnctdnss_ == 1
            [C] = connectedness(Phi(1 : ny*lags, 1 : ny),Sigma,nethor);
        end
        CnndtnssIndex(d,1)         = C.Index;
        CnndtnssFromAlltoUnit(:,d) = C.FromAllToUnit;
        CnndtnssFromUnitToAll(:,d) = C.FromUnitToAll;
        Ctheta(:,:,d)              = C.theta;
    end
    
    % EXOGENOUS BLOCK: restore ny = # of endogenous variables.
    if exogenous_block == 1
        ny = ny0;
    end
    if waitbar_yes, waitbar(d/K, wb); end
    dd = 0; % reset 
end
if waitbar_yes, close(wb); end

%==========================================================================
% [8] CLASSICAL AND PENALIZED ESTIMATES
%   OLS (stored as BVAR.*_ols) and Ridge / Lasso / Elastic-Net (BVAR.Ridge,
%   BVAR.Lasso, BVAR.ElasticNet), with IRFs and information criteria.
%==========================================================================

% --- OLS (classical inference) ---
BVAR.Phi_ols    = varols.B;
BVAR.e_ols      = varols.u;
BVAR.fe_ols     = fe_(fhor,BVAR.e_ols,BVAR.Phi_ols(1 : ny*lags, 1 : ny));
BVAR.Sigma_ols  = 1/(nobs-nk)*varols.u'*varols.u;
% the model with the lowest IC is preferred
[BVAR.InfoCrit.AIC, BVAR.InfoCrit.HQIC, BVAR.InfoCrit.BIC] = IC(BVAR.Sigma_ols, BVAR.e_ols, nobs, nk);
% whiten the reduced form errors
BVAR.Sigma_ols_chol = chol(BVAR.Sigma_ols,'lower');
BVAR.iSigChol      = inv(BVAR.Sigma_ols_chol);
BVAR.white_e_ols   = varols.u * BVAR.iSigChol';

% OLS irf
% with cholesky
BVAR.ir_ols      = iresponse(BVAR.Phi_ols(1 : ny*lags, 1 : ny),BVAR.Sigma_ols,hor,eye(ny));
% with long run
if long_run_irf == 1
    [irlr,Qlr]              = iresponse_longrun(BVAR.Phi_ols(1 : ny*lags, 1 : ny),BVAR.Sigma_ols,hor,lags);
    BVAR.irlr_ols           = irlr;
    BVAR.Qlr_ols(:,:)       = Qlr;
end
% set identified IRF
if set_irf > 0  
    wb = waitbar(0, ['Generating rotations for set-identification - OLS Estimator']);    
    % with sign restrictions
    if signs_irf == 1
        for d1 = 1 : set_irf
            [BVAR.irsign_ols(:,:,:,d1),BVAR.Omega_ols(:,:,d1)] = ...
                iresponse_sign(BVAR.Phi_ols(1 : ny*lags, 1 : ny),BVAR.Sigma_ols,hor,signs);
            waitbar(d1/set_irf, wb);
        end
    end
    % with narrative and sign restrictions
    if narrative_signs_irf == 1
        for d1 = 1 : set_irf
            [BVAR.irnarrsign_ols(:,:,:,d1),BVAR.Omegan_ols(:,:,d1)] = ...
                iresponse_sign_narrative(BVAR.e_ols,BVAR.Phi_ols(1 : ny*lags, 1 : ny),BVAR.Sigma_ols,hor,signs,narrative);
            waitbar(d1/set_irf, wb);
        end
    end
    % with zeros and sign restrictions
    if zeros_signs_irf == 1         
        for d1 = 1 : set_irf
            [BVAR.irzerosign_ols(:,:,:,d1),BVAR.Omegaz_ols(:,:,d1)] = ...
                iresponse_zeros_signs(BVAR.Phi_ols,BVAR.Sigma_ols,hor,lags,var_pos,f,sr);
            waitbar(d1/set_irf, wb);
        end
    end
    % with options.set
    if setrestr_irf == 1
        optS       = setrestr_opt;
        optS.draws = set_irf;
        [BVAR.irfset_ols,BVAR.Omegaset_ols,BVAR.irfset_weights_ols] = ...
            iresponse_set(BVAR.e_ols,BVAR.Phi_ols(1 : ny*lags, 1 : ny),BVAR.Sigma_ols,hor,setrestr,[],optS);
    end
    close(wb)
end
% proxy
if proxy_irf == 1
    inols.res               = BVAR.e_ols;
    inols.Phi               = BVAR.Phi_ols;
    inols.Sigma             = BVAR.Sigma_ols;
    inols.compute_F_stat    = 1;
    tmp_                    = iresponse_proxy(inols);
    BVAR.irproxy_ols(:,:,1) = tmp_.irs';
    BVAR.proxy.F_m          = tmp_.F_m;
    BVAR.proxy.F_m_rob      = tmp_.F_m_rob;
    BVAR.proxy.R2adj_m      = tmp_.R2adj_m;
    BVAR.proxy.data         = options.proxy;
    clear tmp_
end
% with heteroskedasticity
if heterosked_irf == 1
    [BVAR.irheterosked_ols,BVAR.Omegah_ols] = ...
        iresponse_heterosked(BVAR.Phi_ols(1 : ny*lags, 1 : ny),BVAR.e_ols,hor,heterosked_regimes);    
end
% test the normality of the ols VAR residuals (matlab stat toolbox needed)
if  exist('kstest') ==2
    for gg = 1 : ny
        [H,Pv] = kstest(BVAR.white_e_ols(:,gg));
        BVAR.HP(gg,:) = [H,Pv];
        %      H = 0 => Do not reject the null hypothesis at the 5% significance
        %      level. 
    end
else
    BVAR.HP = [];
end
% --- Penalized approaches (regularization): estimates, IRFs, connectedness ---
penalizationstrn = {'Ridge','Lasso','ElasticNet'};
% loop across penalization approaches
% 'name' is both the BVAR.(name) output field and the selector of the
% matching var{Ridge,Lasso,ElasticNet} struct built in [3].
for pp = 1 : 3
    name   = penalizationstrn{pp};
    active = false;
    switch name
        case 'Ridge'
            if Ridge_ == 1, active = true; varPP = varRidge; end
        case 'Lasso'
            if Lasso_ == 1, active = true; varPP = varLasso; end
        case 'ElasticNet'
            if ElasticNet_ == 1, active = true; varPP = varElasticNet; end
    end
    Phi_ = []; Sigma_ = []; u_ = []; Omega_(:,:,1) = eye(ny);  InfoCrit_ = [];
    if active
        % AR param
        BVAR.(name).Phi = varPP.B;
        Phi_            = BVAR.(name).Phi;
        % Error term
        BVAR.(name).e   = varPP.u;
        u_              = BVAR.(name).e;
        % Covariance Matrix
        BVAR.(name).Sigma = 1/(nobs-nk) * varPP.u' * varPP.u;
        Sigma_             = BVAR.(name).Sigma;
        % Recursive IRFs
        BVAR.(name).ir  = iresponse(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,eye(ny));
        % info crit
        [InfoCrit_.AIC, InfoCrit_.HQIC, InfoCrit_.BIC] ...
            = IC(Sigma_, u_, nobs, nk);
        BVAR.(name).InfoCrit = InfoCrit_;
        % IRFs with different identification schemes
        % point identification: LR
        if long_run_irf == 1
            [BVAR.(name).irlr,Omega_] = iresponse_longrun(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,lags);
        end
        % point identification: with heteroskedasticity
        if heterosked_irf == 1
            [BVAR.(name).irheterosked,Omega_] = iresponse_heterosked(Phi_(1 : ny*lags, 1 : ny),u_,hor,heterosked_regimes);
        end
        % set identification: signs
        if set_irf > 0
            wb = waitbar(0, ['Generating rotations for set-identification - ' name ' Estimator']);
        end
        if signs_irf == 1
            for d1 = 1 : set_irf
                [BVAR.(name).irsign_draws(:,:,:,d1), Omega_(:,:,d1)] = ...
                    iresponse_sign(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,signs);
                waitbar(d1/set_irf, wb);
            end
        end
        % set identification: narrative and sign restrictions
        if narrative_signs_irf == 1
            for d1 = 1 : set_irf
                [BVAR.(name).irnarrsign_draws(:,:,:,d1), Omega_(:,:,d1)] = ...
                    iresponse_sign_narrative(u_,Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,signs,narrative);
                waitbar(d1/set_irf, wb);
            end
        end
        % set identification: zeros and sign restrictions
        if zeros_signs_irf == 1
            for d1 = 1 : set_irf
                [BVAR.(name).irzerosign_draws(:,:,:,d1), Omega_(:,:,d1)] = ...
                    iresponse_zeros_signs(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,lags,var_pos,f,sr);
                waitbar(d1/set_irf, wb);
            end
        end
        % set identification: options.set
        if setrestr_irf == 1
            optS       = setrestr_opt;
            optS.draws = set_irf;
            [BVAR.(name).irfset_draws, Omega_, BVAR.(name).irfset_weights] = ...
                iresponse_set(u_,Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,setrestr,[],optS);
        end
        if set_irf>0, close(wb); end
    end
    if cnnctdnss_ == 1 % default identification (Pesaran and Shin)
        BVAR.(name).Connectedness = connectedness(Phi_(1 : ny*lags, 1 : ny),Sigma_,nethor);
    elseif cnnctdnss_ == 2 % customized identification
        Sigma_lower_chol = chol(Sigma_)';
        Omegam           = median(Omega_,3);
        BVAR.(name).Connectedness = connectedness(Phi_(1 : ny*lags, 1 : ny), Sigma_, nethor, Sigma_lower_chol * Omegam);
    end

end
    
%==========================================================================
% [9] STORING THE BAYESIAN RESULTS
%   The last dimension of the *_draws objects indexes the posterior draw.
%==========================================================================

% inference and IRFs
BVAR.Phi_draws    = Phi_draws;          % draws from the autoregressive part
BVAR.Sigma_draws  = Sigma_draws;        % draws from the covariance matrix
BVAR.Sigma_lower_chol_draw = Sigma_lower_chol_draw;
BVAR.ir_draws     = ir_draws;           % draws from the IRF with cholesky
BVAR.irlr_draws   = irlr_draws;         % draws from the IRF with Long Run
BVAR.Qlr_draws    = Qlr_draws;          % Long Run Rotation matrix
BVAR.lags         = lags;               % lags
BVAR.N            = ny;                 % number of variables
BVAR.e_draws      = e_draws;            % residuals
if flag_fe 
    BVAR.fe_draws     = fe_draws;           % insample forecast errors
end
BVAR.posterior    = posterior;
BVAR.prior        = prior;             % priors used
BVAR.logmlike     = log_dnsty;
BVAR.X            = var.X;              % regressors (including dummy if Minnesota)
BVAR.y            = var.y;              % dependent  (including dummy if Minnesota)
BVAR.XX           = XX;                 % regressors (no dummy)
BVAR.yy           = yy;                 % dependent  (no dummy)

% prediction
BVAR.fhor         = fhor;               % forecast horizon
BVAR.hor          = hor;                % IRF horizon
BVAR.forecasts.no_shocks      = squeeze(yhatfut_no_shocks);         % trajectories of forecasts without shocks
BVAR.forecasts.with_shocks    = squeeze(yhatfut_with_shocks);       % trajectories of forecasts with shocks
BVAR.forecasts.conditional    = [];                                 % trajectories of conditional forecasts
BVAR.forecasts.EPScond        = [];                                 % shocks of conditional forecasts
if cfrcst_yes ~= 0
    BVAR.forecasts.conditional    = squeeze(yhatfut_cfrcst);   % trajectories of forecasts
    BVAR.forecasts.EPScond        = EPS;                       % shocks of forecasts
end
BVAR.forecast_data            = forecast_data;
if robust_bayes_ > 0
    BVAR.Kstar_  = Kstar_;
    if robust_bayes_ > 1
        BVAR.Sstar_  = Sstar_;
    end
    BVAR.hom_shrinkage = K_shrinkage/(nobs+K_shrinkage);
end
%
BVAR.varnames     = varnames;
BVAR.ndraws       = K;

if signs_irf == 1 && narrative_signs_irf == 0
    BVAR.irsign_draws = irsign_draws;
    BVAR.Omegas       = Omega_draws;
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irsign_lower_draws = irsign_lower_draws;
        BVAR.irsign_upper_draws = irsign_upper_draws;
        BVAR.irsign_robust_credible_bands = compute_rob_cred_bands_(irsign_lower_draws, irsign_upper_draws, ny, opt_GiacomoniKitagawa);
    end
else
    BVAR.irsign_draws = [];
    BVAR.Omegas       = [];
end
if narrative_signs_irf == 1
    BVAR.irnarrsign_draws = irnarrsign_draws;
    BVAR.Omegan           = Omegan_draws;
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irnarrsign_lower_draws = irnarrsign_lower_draws;
        BVAR.irnarrsign_upper_draws = irnarrsign_upper_draws;
        BVAR.irnarrsign_robust_credible_bands = compute_rob_cred_bands_(irnarrsign_lower_draws, irnarrsign_upper_draws, ny, opt_GiacomoniKitagawa);
    end    
else
    BVAR.irnarrsign_draws = [];
    BVAR.Omegan           = [];
end
if zeros_signs_irf == 1
    BVAR.irzerosign_draws   = irzerosign_draws;
    BVAR.Omegaz             = Omegaz_draws;
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irzerosign_lower_draws = irzerosign_lower_draws;
        BVAR.irzerosign_upper_draws = irzerosign_upper_draws;
        BVAR.irzerosign_robust_credible_bands = compute_rob_cred_bands_(irzerosign_lower_draws, irzerosign_upper_draws, ny, opt_GiacomoniKitagawa);
    end
else
    BVAR.irzerosign_draws = [];
    BVAR.Omegaz           = [];
end
if setrestr_irf == 1
    BVAR.irfset_draws       = irfset_draws;
    BVAR.Omegaset           = Omegaset_draws;
    BVAR.irfset_weights     = irfset_weights;   % importance weights (0 = no rotation found)
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irfset_lower_draws = irfset_lower_draws;
        BVAR.irfset_upper_draws = irfset_upper_draws;
        BVAR.irfset_robust_credible_bands = compute_rob_cred_bands_(irfset_lower_draws, irfset_upper_draws, ny, opt_GiacomoniKitagawa);
    end
else
    BVAR.irfset_draws   = [];
    BVAR.Omegaset       = [];
    BVAR.irfset_weights = [];
end
if heterosked_irf == 1
    BVAR.irheterosked_draws   = irheterosked_draws;
    BVAR.Omegah_draws         = Omegah_draws;
    BVAR.Omegah               = BVAR.Omegah_draws; 
else
    BVAR.irheterosked_draws   = [];
    BVAR.Omegah_draws         = [];
    BVAR.Omegah               = [];
end
if proxy_irf == 1
    BVAR.irproxy_draws = irproxy_draws;
    BVAR.Omegap        = Omegap_draws;
else
    BVAR.irproxy_draws= [];
    BVAR.Omegap = [];
end
if hmoments_signs_irf == 1
    BVAR.irhmomsign_draws = irhmomsign_draws;
    BVAR.Omegam           = Omegam_draws;
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irhmomsign_lower_draws = irhmomsign_lower_draws;
        BVAR.irhmomsign_upper_draws = irhmomsign_upper_draws;
        BVAR.irhmomsign_robust_credible_bands = compute_rob_cred_bands_(irhmomsign_lower_draws, irhmomsign_upper_draws, ny, opt_GiacomoniKitagawa);
    end

else
    BVAR.irhmomsign_draws = [];
    BVAR.Omegam           = [];
end
if hmoments_eig_irf == 1
    BVAR.irhmomeig_draws = irhmomeig_draws;
    BVAR.Omegae           = Omegae_draws;
else
    BVAR.irhmomeig_draws = [];
    BVAR.Omegae          = [];
end
if nexogenous > 0
    BVAR.irx_draws = irx_draws;
else
    BVAR.irx_draws = [];
end

% missing observations
BVAR.data         = y;                  % raw data
if mixed_freq_on
    BVAR.yfill = yfill;
    BVAR.yfilt = yfilt;
    BVAR.yinterpol = yinterpol;
    BVAR.logL  = logL;
    if isfield(options,'mf_varindex')== 1
        KFoptions.mf_varindex = options.mf_varindex;
    end
    BVAR.KFoptions        =  KFoptions;
end

% bvar connectedness
if cnnctdnss_
    BVAR.Connectedness.Index         = CnndtnssIndex;
    BVAR.Connectedness.FromAlltoUnit = CnndtnssFromAlltoUnit;
    BVAR.Connectedness.FromUnitToAll = CnndtnssFromUnitToAll;
    BVAR.Connectedness.theta         = Ctheta;
end


%==========================================================================
% [10] NESTED FUNCTION posterior_
%==========================================================================
    function [posterior,var] = posterior_(y)
        
        % POSTERIOR_ Posterior moments of the VAR given data y.
        %   posterior.PhiHat, .XXi, .S, .df   Normal-Inverse-Wishart moments
        %   var                               rfvar3 output (with prior dummies)
        %   Called (i) once in [5]; (ii) for the z-block of options.exogenous_block;
        %   (iii) at each draw with mixed frequency, on the Kalman-smoothed data.
        %   NESTED FUNCTION: it shares the workspace of bvar_. It READS idx, nx,
        %   ny, lags, timetrend, nexogenous, exogenous, pandemic_on, panDum, ydum,
        %   xdum, pbreaks, lambda, mu, plr, ww, heterosked, nunits, dummy, flat,
        %   prior, exogenous_block, nz, XX2, and it OVERWRITES ydata, T and xdata
        %   of bvar_ (hence ydata0 in [5]).
        %===========================================
        ydata = y(idx, :, :);
        %===========================================
        T     =  size(ydata, 1);
        xdata = ones(T,nx);
        if timetrend ==1
            xdata = [xdata [1-lags:T-lags]'];
        end
        if nexogenous >0
            xdata = [xdata exogenous(idx,:)];
        end

        if pandemic_on == 1
           xdata = [xdata, panDum(idx, :)];
        end


        % EXOGENOUS BLOCK: y-block regression -> append the z lags. When called
        % for the z-block, the flag is temporarily 0 so z is regressed on its own
        % lags and the constant only.
        if exogenous_block == 1
            xdata = [xdata,  [nan(lags,lags*nz); XX2(:,1:lags*nz)]];
        end
        % posterior density
        if nunits == 1
            if heterosked == 0
                var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu, [], plr);
            else
                var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu, ww, plr);
            end
        else  % pooled units
            var.y = []; var.X = [];
            for nunt = 1 : nunits
                tmp_var = rfvar3(ydata(:,:,nunt), lags, xdata, [T; T], 0, 0);
                var.X = [var.X; tmp_var.X];
                var.y = [var.y; tmp_var.y];
            end
            % Compute OLS regression and residuals of the pooled estimator
            var = fast_ols_(var.y,var.X);
        end
        Tu = size(var.u, 1);
                
        if dummy == 1
            %********************************************************
            % Minnesota Prior
            % prior.df/S/XXi/PhiHat are computed once before the MCMC
            % loop (in the prior-specification block above posterior_).
            % They are fixed across draws: the presample is assumed to
            % be complete (no missing values), so the Kalman-smoothed
            % data passed as y does not affect the prior density.
            %********************************************************
            posterior.df    = Tu - ny*lags - nx - flat*(ny+1);
            posterior.S     = var.u' * var.u;
            posterior.XXi   = var.xxi;
            posterior.PhiHat = var.B;
            
        elseif dummy == 2
            %********************************************************
            % Conjugate Prior
            %********************************************************
            Ai              = inv(prior.Phi.cov);
            posterior.df    = Tu + prior.Sigma.df;
            posterior.XXi   = inv(var.X'*var.X + Ai);
            posterior.PhiHat = posterior.XXi * (var.X' * var.y + Ai * prior.Phi.mean);
            %posterior.S     = var.u' * var.u + prior.Sigma.scale ;
            %posterior.S     = var.u' * var.u + prior.Sigma.scale  + (var.B - prior.Phi.mean)' * Ai * (var.B - prior.Phi.mean);
            %posterior.S     = var.u' * var.u + prior.Sigma.scale  + ...
            %    (posterior.PhiHat - prior.Phi.mean)' * Ai * (var.B - prior.Phi.mean);
            % check 
            posterior.S = var.u' * var.u + prior.Sigma.scale + ...
                prior.Phi.mean' * Ai * prior.Phi.mean + ...
                var.B' * (var.X'*var.X) * var.B ...
                - posterior.PhiHat' * (var.X'*var.X + Ai) * posterior.PhiHat;
            
        elseif dummy == 0
            %********************************************************
            % Flat Jeffrey Prior
            %********************************************************
            posterior.df    = Tu - ny*lags - nx + flat*(ny+1) - nexogenous;
            posterior.S     = var.u' * var.u;
            posterior.XXi   = var.xxi;
            posterior.PhiHat = var.B;
        

        elseif dummy == 4
            %********************************************************
            % Pandemic Minnesota Prior
            %********************************************************
            posterior.df     = Tu - ny*lags - nx - flat*(ny+1) - pandemic_h;
            posterior.S      = var.u' * var.u;
            posterior.XXi    = var.xxi;
            posterior.PhiHat = var.B;

        elseif dummy == 3
            %********************************************************
            % PLR prior (Minnesota short-run + PLR long-run)
            % prior.df/S/XXi/PhiHat are computed once before the MCMC
            % loop (in the prior-specification block above posterior_).
            % They are fixed across draws: the presample is assumed to
            % be complete (no missing values), so the Kalman-smoothed
            % data passed as y does not affect the prior density.
            %********************************************************
            posterior.df    = Tu - ny*lags - nx - flat*(ny+1);
            posterior.S     = var.u' * var.u;
            posterior.XXi   = var.xxi;
            posterior.PhiHat = var.B;
        end
         
    end

% end of bvar_ (posterior_ above is nested in it)
end
