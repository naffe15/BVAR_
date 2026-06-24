function opt = parse_bvar_options(y, lags, options)
%PARSE_BVAR_OPTIONS Default settings, options parsing, and mixed-frequency
% interpolation for bvar_.
%
% opt = parse_bvar_options(y, lags, options)
% opt = parse_bvar_options(y, lags)              % no options given
%
% This is a near-verbatim relocation of what used to be the first ~640
% lines of bvar_.m (everything between the nargin/lags validity checks
% and the "Consistency Checks" section). It existed there because that's
% where the data first arrives, not because the logic itself depends on
% anything else in bvar_ -- so it was safe to move as a unit. Nothing in
% the ~1300 lines that follow in bvar_.m had to change: every variable
% this function used to set as a bare local is returned here as a field
% of OPT instead, and bvar_.m unpacks them back into the same names right
% after calling this.
%
% A few things this function does that are worth knowing if you're
% tracing a bug back to here:
%   - It can MUTATE the data: if y contains NaNs, it's treated as mixed-
%     frequency/irregularly sampled and interpolated in place (opt.y is
%     the interpolated series; opt.yoriginal is the unmodified input).
%   - 'prior' (singular) and 'priors' (plural) are different variables.
%     'priors' is just a name tag ('Minnesota'/'Conjugate'/'N/A'). 'prior'
%     only gets populated here under the Conjugate-prior branch
%     (prior.Phi.mean/cov, prior.Sigma.scale/df) and isn't read again
%     until deep inside bvar_'s posterior computation -- easy to miss if
%     you're scanning for where it's used.
%   - The original code called a tiny nested function 'priors_()' (just
%     `priors.name = 'N/A'`) defined inside bvar_.m. Nested functions
%     can't be called from a separate file, so that single line is
%     inlined directly below instead -- behaviorally identical.
%
% See bvar_.m for the full list of supported options.* fields.

% number of observable variables
ny                  = size(y, 2);
% number of units (for panels)
nunits              = size(y, 3);

%********************************************************
%* DEFAULT SETTINGS
%********************************************************
% Control random number generator
if isOctave == 0
    isMatlab = 1;
    rng('default');
    rng(999);
else
    isMatlab = 0;
    % pkg load optim
    randn('state',999);
    rand('state',999);
end


% Default Settings (they can all be changed in 'options' see below)
K                   = 5000;         % number of draws from the posterior
hor                 = 24;           % horizon for the IRF
fhor                = 12;           % horizon for the forecasts
nethor              = 12;           % network horizon/ splioover and connectedness
firstobs            = lags+1;       % first observation
presample           = 0;            % using a presample for setting the hyper-parameter of the Minnesosta prior
noconstant          = 0;            % when 0, includes a constatn in the VAR
timetrend           = 0;            % when 1, includes a time trend in the VAR
minn_prior_tau      = 3;            % Minnesota prior Hyper-Param: Overall Tightness
minn_prior_decay    = 0.5;          % Minnesota prior Hyper-Param: Tighness on lags>1
minn_prior_lambda   = 5;            % Minnesota prior Hyper-Param: Sum-of-Coefficient
minn_prior_mu       = 2;            % Minnesota prior Hyper-Param: Co-Persistence
minn_prior_omega    = 2;            % Minnesota prior Hyper-Param: Shocks Variance
long_run_irf        = 0;            % when 0, it does not compute long run IRF
irf_1STD            = 1;            % when 1, IRF are computed as 1SD increase. Else, IRF are compued as unitary increase in the shock
cfrcst_yes          = 0;            % no conditional forecast unless defined in options
non_explosive_      = 0;            %
heterosked          = 0;

signs_irf           = 0;
narrative_signs_irf = 0;
zeros_signs_irf     = 0;
proxy_irf           = 0;
heterosked_irf      = 0;
hmoments_signs_irf  = 0;
hmoments_eig_irf    = 0;
noprint             = 0;
nexogenous          = 0;
exogenous           = [];
cnnctdnss_          = 0;
Ridge_              = 0;
Lasso_              = 0;
ElasticNet_         = 0;
set_irf             = 0;
robust_bayes_       = 0;
robust_credible_regions_  = 0;
exogenous_block = 0;
nz              = 0;

% for mixed frequecy / irregurerly sampled data.
% Interpolate the missing values of each times series.
mixed_freq_on = 0;
if any(any(isnan(y))) %== 1
        mixed_freq_on = 1;
    index_nan_var = find(sum(isnan(y),1) ~= 0);
    yoriginal     = y;
    % interporlate the variables that have nans
    T = 1:1:length(y);
    for kk  = 1 : length(index_nan_var)
        v                       = y(isfinite(y(:,index_nan_var(kk))),index_nan_var(kk));
        x                       = find(isfinite(y(:,index_nan_var(kk))));
        if isMatlab == 1
            y(:,index_nan_var(kk))  = interp1(x,v,T','spline');
        else
            y(:,index_nan_var(kk))  = spline(x,v,T');
        end
        yinterpol               = y;
    end
end
if mixed_freq_on == 1 && nargin < 3
    warning('You did not specified the aggregation of the mixed freq. Variables will be treated as stocks.');
    index = zeros(size(y,2),1);
end
if mixed_freq_on == 1 && nunits > 1
    error('The toolbox does not estimate a pooled VAR with missing observations.');
end

% Priors declaration: default Jeffrey prior
dummy    = 0;
flat     = 1;
priors.name  = 'N/A';   % was priors_() -- see file header note

% declaring the names for the observable variables
for v = 1 : ny
    varnames{v} = ['Var' num2str(v)];
end


%********************************************************
%* CUSTOMIZED SETTINGS
%********************************************************
if nargin > 2
    if isfield(options,'vnames')==1
        varnames = options.vnames;
    end
    %======================================================================
    % Inference options
    %======================================================================
    if isfield(options,'K')==1
        K = options.K;
    end
    if isfield(options,'firstobs')==1
        firstobs = options.firstobs;
        if firstobs < lags + 1
            error('firstobs need to be larger than lags +1')
        end
    end
    if isfield(options,'presample')==1
        presample = options.presample;
    end
    if isfield(options,'noconstant')==1
        noconstant = options.noconstant;
    end
    if isfield(options,'timetrend')==1
        timetrend  = options.timetrend;
        noconstant = 0;
    end
    if isfield(options,'non_explosive_')==1
        non_explosive_  = options.non_explosive_;
    end
    if isfield(options,'heterosked_weights')==1
        ww  = options.heterosked_weights;
        heterosked = 1;
    end
    % compute the second and fourth moments robust to error distribution misspecification
    if isfield(options,'robust_bayes')==1
        robust_bayes_  = options.robust_bayes;
        % shrinkage towards a normal distribution K_shrinkage = infty
        if isfield(options,'K_shrinkage') == 1
            K_shrinkage = options.K_shrinkage;
        else
            K_shrinkage = 1 * size(y, 1);
        end
    end
    %======================================================================
    % Exogenous Variables options
    %======================================================================
    if isfield(options,'exogenous')==1 || isfield(options,'controls')==1
        if isfield(options,'controls')==1
            exogenous = options.controls;
        else
            exogenous = options.exogenous;
        end
        nexogenous = size(exogenous,2);
        if size(exogenous,1) ~= size(y,1)+ fhor && size(exogenous,1) ~= size(y,1)
            error('Size Mismatch between endogenous and exogenos variables; exo must be either T or T+fhor');
        end
        if any(isnan(exogenous(lags+1:end,:)))
            error('Exogenous variables cannnot be ''nan'' from lags+1 onward.');
        end
    end
    %======================================================================
    % Exogenous Block options
    %======================================================================
    if isfield(options,'exogenous_block')==1
        exogenous_block = 1;
        exogenous = options.exogenous_block;

        nz = size(exogenous,2);
        if size(exogenous,1) ~= size(y,1)+ fhor && size(exogenous,1) ~= size(y,1)
            error('Size Mismatch between endogenous and exogenos variables; exo must be either T or T+fhor');
        end
        if any(isnan(exogenous(lags+1:end,:)))
            error('Exogenous variables cannot be ''nan'' from lags+1 onward.');
        end
    end
    %======================================================================
    % Minnesota prior options
    %======================================================================
    if (isfield(options,'priors')==1 && strcmp(options.priors.name,'Minnesota')==1) || (isfield(options,'priors')==1 && strcmp(options.priors.name,'minnesota')==1) || ...
       (isfield(options,'prior')==1 && strcmp(options.prior.name,'Minnesota')==1) || (isfield(options,'prior')==1 && strcmp(options.prior.name,'minnesota')==1)
        %  MINNESOTA PRIOR
        dummy = 1;
        flat  = 0;
        timetrend = 0;
        priors.name= 'Minnesota';
    end
    if isfield(options,'minn_prior_tau')==1 || isfield(options,'bvar_prior_tau')==1
        dummy = 1;
        flat  = 0;
        priors.name= 'Minnesota';
        %  MINNESOTA PRIOR: tightness
        if isfield(options,'bvar_prior_mu')==1
            minn_prior_tau = options.bvar_prior_tau;
        else
            minn_prior_tau = options.minn_prior_tau;
        end
    end
    if isfield(options,'minn_prior_decay')==1 || isfield(options,'bvar_prior_decay')==1
        dummy = 1;
        flat  = 0;
        priors.name= 'Minnesota';
        %  MINNESOTA PRIOR: decay
        if isfield(options,'bvar_prior_mu')==1
            minn_prior_decay = options.bvar_prior_decay;
        else
            minn_prior_decay = options.minn_prior_decay;
        end
    end
    if isfield(options,'minn_prior_lambda')==1 || isfield(options,'bvar_prior_lambda')==1
        dummy = 1;
        flat  = 0;
        priors.name= 'Minnesota';
        %  MINNESOTA PRIOR: sum-of-coeff
        if isfield(options,'bvar_prior_mu')==1
            minn_prior_lambda = options.bvar_prior_lambda;
        else
            minn_prior_lambda = options.minn_prior_lambda;
        end
    end
    if isfield(options,'minn_prior_mu')==1 || isfield(options,'bvar_prior_mu')==1
        dummy = 1;
        flat  = 0;
        priors.name= 'Minnesota';
        %  MINNESOTA PRIOR: co-persistence
        if isfield(options,'bvar_prior_mu')==1
            minn_prior_mu = options.bvar_prior_mu;
        else
            minn_prior_mu = options.minn_prior_mu;
        end
    end
    if isfield(options,'minn_prior_omega')==1 || isfield(options,'bvar_prior_omega')==1
        dummy = 1;
        flat  = 0;
        priors.name= 'Minnesota';
        %  MINNESOTA PRIOR: variance
        if isfield(options,'bvar_prior_omega')==1
            minn_prior_omega = options.bvar_prior_omega;
        else
            minn_prior_omega = options.minn_prior_omega;
        end
    end
    if isfield(options,'max_minn_hyper')==1 && options.max_minn_hyper ==1 && mixed_freq_on ==0
        % maximize the hyper parameters of the minnesota prior
        hyperpara(1) = minn_prior_tau;
        hyperpara(2) = minn_prior_decay;
        hyperpara(3) = minn_prior_lambda;
        hyperpara(4) = minn_prior_mu;
        hyperpara(5) = minn_prior_omega;

        dummy = 1;
        flat  = 0;
        priors.name= 'Minnesota';
        try
            [postmode,lm,~] = bvar_max_hyper(hyperpara,y,lags,options);

            minn_prior_tau      = postmode(1);
            minn_prior_decay    = postmode(2);
            minn_prior_lambda   = postmode(3);
            minn_prior_mu       = postmode(4);
            minn_prior_omega    = postmode(5);
            disp('                       ')
            disp('Maximization Successful: I will use the mode values.')

        catch
            warning('Maximization NOT Successful')
            disp('Using hyper parameter default values')
        end
    end
    %======================================================================
    % Conjugate/Hierachical MN-IW prior options
    %======================================================================
    if (isfield(options,'priors')==1 && strcmp(options.priors.name,'Conjugate')==1) || (isfield(options,'priors')==1 && strcmp(options.priors.name,'conjugate')==1) || ...
       (isfield(options,'prior')==1 && strcmp(options.prior.name,'Conjugate')==1) || (isfield(options,'prior')==1 && strcmp(options.prior.name,'conjugate')==1)

        if isfield(options,'prior')==1
            options.priors = options.prior;
        end
        if dummy == 1
            warning('You have set both the Conjugate and Minnesota (perhaps via options.minn_prior_XX)');
            warning('I will consider the Conjugate prior only');
        end
        dummy = 2;
        %  warning('The Conjugate prior is still under construction ... ');
        flat  = 0;
        priors.name= 'Conjugate';
        % Priors for the AR parameters
        if isfield(options.priors,'Phi') == 1
            % mean
            if isfield(options.priors.Phi,'mean') == 1
                prior.Phi.mean  = options.priors.Phi.mean;
                if max(size(prior.Phi.mean) ~= [ny*lags+(1-noconstant)+timetrend+nexogenous   ny]) ~= 0
                    error('Size mismatch')
                end
            else
                warning(['You did not provide a prior mean for the AR coeff. ' ...
                    'Assume zeros everywhere.'])
                prior.Phi.mean  = zeros(ny*lags+(1-noconstant)+timetrend+nexogenous , ny);
            end
            % variance
            if isfield(options.priors.Phi,'cov') == 1
                prior.Phi.cov   = options.priors.Phi.cov;
%                 if length(prior.Phi.cov) ~= (ny*lags+(1-noconstant) +timetrend ) * ny
                if length(prior.Phi.cov) ~= (ny*lags+(1-noconstant) +timetrend+nexogenous ) || size(prior.Phi.cov,1 )~=size(prior.Phi.cov,2)
                    error('Size mismatch: Covariance Phi should be square, e.g. size(Phi.mean,1)x size(Phi.mean,1)')
                end
            else
                warning(['You did not provide a Covariance for the AR coeff. ' ...
                    'Assume 10 times Identity Matrix.'])
%                 prior.Phi.cov  = 10 * eye((ny*lags+(1-noconstant) + timetrend) * ny);
                prior.Phi.cov  = 10 * eye((ny*lags+(1-noconstant) + timetrend+nexogenous));
            end
        else
            warning(['You did not provide prior mean and covariance for the AR coeff ' ...
                'Assume zeros everywhere with covariance 10 times Identity Matrix.'])
%             prior.Phi.cov   = 10 * eye((ny*lags+(1-noconstant) +timetrend ) * ny);
            prior.Phi.cov   = 10 * eye((ny*lags+(1-noconstant) + timetrend+nexogenous ));
            prior.Phi.mean  = zeros(ny*lags+(1-noconstant) + timetrend+nexogenous ,  ny);
        end
        % Priors for the Residual Covariance
        if isfield(options.priors,'Sigma') == 1
            % scale
            if isfield(options.priors.Sigma,'scale') == 1
                prior.Sigma.scale = options.priors.Sigma.scale;
                if size(prior.Sigma.scale) ~= [ny ny]
                    error('Size mismatch')
                end
            else
                warning(['You did not provide a prior scale for the Residual Covariance. ' ...
                    'Assume identity matrix.'])
                prior.Sigma.scale = eye(ny);
            end
            % degrees of freedom
            if isfield(options.priors.Sigma,'df') == 1
                prior.Sigma.df = options.priors.Sigma.df;
                if length(prior.Sigma.df) ~= 1
                    error('Size mismatch')
                end
                if prior.Sigma.df/2 <= ny-1
                    error('Too few degrees of freedom - Increase prior df')
                end
            else
                warning(['You did not provide the degrees of freedom for the Residual Covariance. ' ...
                    'Assume N+1 degrees of freedom.'])
                prior.Sigma.df = ny + nexogenous + timetrend + 1;
                while prior.Sigma.df/2 <= ny-1 % too few df
                    prior.Sigma.df =  prior.Sigma.df +1;
                end
            end
        else
            warning(['You did not provide prior scale and degrees of freedom for the Residual Covariance. ' ...
                'Assume an identity matrix matrix with N+1 degrees of freedom.'])
            prior.Sigma.scale = eye(ny);
            prior.Sigma.df    = ny + nexogenous + timetrend + 1;
            while prior.Sigma.df/2 <= ny-1 % too few df
                prior.Sigma.df =  prior.Sigma.df +1;
            end
        end
    end
    %======================================================================
    % IRF options
    %======================================================================
    if isfield(options,'hor') ==1
        hor = options.hor;
    end
    if isfield(options,'set_irf') ==1
        % set_irf is the # of rotations for set identified systems in wihch
        % Phi and Sigma are point estimates (e.g. OLS, Ridge, Lasso or
        % Elastic Nets); default 0
        set_irf = options.set_irf;
    end
    if isfield(options,'robust_credible_regions') ==1
        % computing the robust bayesian credible regions of
        % Giacomini-Kitagawa ECMA 2021
        robust_credible_regions_ = 1;
        if isfield(options.robust_credible_regions,'KG') == 1 && options.robust_credible_regions.KG == 0 % to disactivate
            robust_credible_regions_ = 0;
        end
        L                        = 1000;
        opt_GiacomoniKitagawa.aalpha    = 0.68; % Credibility level
        opt_GiacomoniKitagawa.gridLength = 1000;
        if isfield(options.robust_credible_regions,'L') == 1
            L = options.robust_credible_regions.L;
        end
        if isfield(options.robust_credible_regions,'aalpha') == 1
            opt_GiacomoniKitagawa.aalpha = options.robust_credible_regions.aalpha;
        end
        if isfield(options.robust_credible_regions,'gridLength') == 1
            opt_GiacomoniKitagawa.gridLength = options.robust_credible_regions.gridLength;
        end
    end
    if isfield(options,'long_run_irf')==1
        % Activating Long run IRF
        long_run_irf = options.long_run_irf;
    end
    if isfield(options,'heterosked_regimes')==1
        % Activating identification via heteroskedasticity
        if length(options.heterosked_regimes) ~= size(y,1)
            error('heterosked_regimes must have the same time dimension of y')
        end
        heterosked_irf     = 1;
        heterosked_regimes = options.heterosked_regimes(lags+1:end);
        if any(heterosked_regimes>1) || any(heterosked_regimes<0)
            error('heterosked_regimes must contain zero (first regime) and one (second regime)')
        end
    end
    if isfield(options,'signs')==1
        % Activating IRF with sign restrictions (mulitple horizons allowed)
        signs_irf       = 1;
        signs           = options.signs;
        if iscellstr(signs) == 0
            error(['options.signs should be a cell array. Each cell must contain a string with the format'...
                '\n''y(a,b,c)<0'' or ''y(a,b,c)>0'' where a, b and c are integers.',...
                '\na = index of the variable',...
                '\nb = horizon',...
                '\nc = index of the shock'],class(zeros))
        end
    end
    if isfield(options,'narrative')==1
        if signs_irf  == 0
            warning('You did not provide any sign restrictions.')
            signs{1} = 'isempty(y(1,1,1))==0';
        end
        if signs_irf  == 1
            signs_irf       = 0;  % disactivating signs
        end
        narrative_signs_irf = 1;
        narrative           = options.narrative ;
    end
    if isfield(options,'zeros_signs')==1
        if signs_irf  == 1
            signs_irf       = 0;  % disactivating signs
        end
        if narrative_signs_irf  == 1
            narrative_signs_irf       = 0;  % disactivating narrative
        end
        % Activating IRF with zeros and sign restrictions (mulitple horizons NOT allowed)
        zeros_signs_irf  = 1;
        zeros_signs      = options.zeros_signs;
        if iscellstr(zeros_signs) == 0
            error(['options.zeros_signs should be a cell array.'...
                '\nEach cell must contain a string with the following format'...
                '\nFor sign restrictions ''y(a,b)=1'' or ''y(a,b)=-1'',',...
                '\nFor short run zero restriction ''ys(a,b)=0'',',...
                '\nFor long run restriction ''yr(a,1,b)=0'' where a and b are integers.',...
                '\na = index of the variable',...
                '\nb = index of the shock'],class(zeros_signs))
        end
        [f,sr] = sign2matrix(zeros_signs,ny+nz);
        if isfield(options,'var_pos')==0
            var_pos = ones(1,ny);
        else
            var_pos = options.var_pos;
        end
    end
    if isfield(options,'proxy')==1
        % Activating IRF with provy
        proxy_irf   = 1;
        in.proxies  = options.proxy;
        in.vars     = y;
        in.p        = lags;
        in.compute_F_stat = 0;
        if isfield(options,'proxy_end') == 1
            in.T_m_end  = options.proxy_end;
        else
            in.T_m_end  = 0;  %if the times series of the instrument ends when VAR data ends
        end
        in.irhor    = hor;
        if isnumeric(in.proxies) == 0
            error(['options.proxy should be a numeric array (nans or inf not allowed)'],class(in.proxies))
        end
        inols =in;
    end
    if isfield(options,'hmoments')==1
        if signs_irf  == 0
            warning('You did not provide any sign restrictions.')
            signs{1} = 'isempty(y(1,1,1))==0';
        end
        if signs_irf  == 1
            signs_irf       = 0;  % disactivating signs
        end
        if robust_bayes_ == 0
            robust_bayes_      = 1;
            K_shrinkage        = 1 * size(y,2);
        end
        hmoments_signs_irf = 1;
        hmoments           = options.hmoments;
        [f]                = hmoments2matrix(hmoments,ny);
    end
    if isfield(options,'hmoments_eig')==1
        moment = options.hmoments_eig;
%         if moment ~= 3 ||  moment ~= 4
%             warning('You can decompose third (3) or fourth (4) moments. I set as default 4')
%             moment = 4;
%         end
        if robust_bayes_ == 0
            robust_bayes_      = 1;
            K_shrinkage        = 1 * size(y,2);
        end
        hmoments_eig_irf = 1;
    end
    if isfield(options,'irf_1STD')==1
        % Activating of unitary IRF, i.e. a unitary increase in the shocks
        % (instead of 1 STD)
        irf_1STD = options.irf_1STD;
    end
    %======================================================================
    % (Un)Conditional Forecasts options
    %======================================================================
    if isfield(options,'fhor')==1
        fhor = options.fhor;
        if fhor < 1
            error('Forecast horizon must be positive')
        end
    end
    if isfield(options,'endo_index')==1
        % Forecast conditional on the path of an endogenous var
        cfrcst_yes      = 1;
        if isfield(options,'endo_path')== 0
            error('You need to provide the p[ath for the endogenou variable')
        end
        % rows forecasts, column variables
        endo_path       = options.endo_path;
        endo_path_index = options.endo_index;
        if length(endo_path_index) ~= size(endo_path)
            error(['Mismatch beween the number of endogenous paths and the number of conditioned variables'...
                '\nE.g. the # of conditioned variables must coincide with the # of column in ''options.endo_path'''],class(endo_path_index));
        end
        if isfield(options,'exo_index')==1
            % Forecast conditional on the path of an endo var using only a
            % subset of shocks. notice that the # of endo and # exo must coincide
            cfrcst_yes  = 2;
            exo_index   = options.exo_index;
            Omega       = eye(ny);
            %             if isfield(options,'Omegaf')==1
            %                 Omegaf = options.Omegaf;
            %             end
            if length(exo_index) ~= length(endo_path_index)
                error('the # of conditioned endogenous and # exogenous shocks used must coincide');
            end
        end
        if nunits > 1
            cfrcst_yes = 0;
            warning('Conditional forecasts are not supported with pooled VARs.')
        end
    end
    %======================================================================
    % Missing Values options
    %======================================================================
    if mixed_freq_on == 1 && isfield(options,'mixed_freq_index')==1
        index = options.mixed_freq_index;
        if length(index) ~= size(y,2)
            error(['You have to specify as many index as observables.'...'
                '\nIf no missing values for var j index(j)=0'...
                '\nIf missing values for var j, and var j is a stock index(j)=0'...
                '\nIf missing values for var j, and var j is a real flow index(j)=2'],class(index))
        end
    elseif mixed_freq_on == 1 && isfield(options,'mf_varindex')== 1
        index = zeros(ny,1);
        index(options.mf_varindex) = 2;
    elseif (mixed_freq_on == 1 && isfield(options,'mixed_freq_index')== 0) || (mixed_freq_on == 1 && isfield(options,'mf_varindex')== 0)
        warning(['You did not specified the aggregation of the mixed freq. Variables will be treated as stocks.']);
        index = zeros(size(y,2),1);
    end
    if isfield(options,'noprint')==1
        noprint = options.noprint;
    end
    %======================================================================
    % Network-Spillover-Connectedness options
    %======================================================================
    if isfield(options,'nethor')==1
        nethor = options.nethor;
    end
    if isfield(options,'connectedness')==1
        cnnctdnss_ = options.connectedness;
        if cnnctdnss_ == 2 && set_irf == 0
           set_irf = 100;
        end
    end
    if isfield(options,'Ridge')==1
        Ridge_     = 1;%options.Ridge;
        if isfield(options.Ridge,'est') == 1
            Ridge_ = options.Ridge.est;
        end
        %cnnctdnss_ = 1;
        if isfield(options.Ridge,'lambda') == 1
            Ridge_lambda = options.Ridge.lambda;
        else
            warning('You did not specify a value for the penalization parameter (options.Ridge.lambda)')
            warning('I am using lambda = 0.02')
            Ridge_lambda = 0.02;
        end
    end
    if isfield(options,'Lasso')==1
        % (matlab stat toolbox needed)
        if exist('lasso') ~= 2
            error('Cannot estimate VAR with Lasso: matlab stat toolbox needed')
        end
        Lasso_     = 1;%options.Lasso;
        if isfield(options.Lasso,'est') == 1
            Lasso_ = options.Lasso.est;
        end
        %cnnctdnss_ = 1;
        if isfield(options.Lasso,'lambda') == 1
            Lasso_lambda = options.Lasso.lambda;
        else
            warning('You did not specify a value for the penalization parameter (options.Lasso.lambda)')
            %warning('Use the largest value of Lambda that gives a nonnull model')
            warning('I am using lambda = 0.05')
            Lasso_lambda = 0.05;
        end
    end
    if isfield(options,'ElasticNet')==1
        % (matlab stat toolbox needed)
        if  exist('lasso') ~= 2
            error('Cannot estimate VAR with ElasticNet: matlab stat toolbox needed')
        end
        ElasticNet_ = 1;%options.ElasticNet;
        if isfield(options.ElasticNet,'est') == 1
            ElasticNet_ = options.ElasticNet.est;
        end
        %cnnctdnss_  = 1;
        if isfield(options.ElasticNet,'lambda') == 1
            ElasticNet_lambda = options.ElasticNet.lambda;
        else
            warning('You did not specify a value for the penalization parameter (options.ElasticNet.lambda)')
            warning('I am using lambda = 0.05')
            ElasticNet_lambda = 0.05;
        end
        if isfield(options.ElasticNet,'alpha') == 1
            ElasticNet_alpha = options.ElasticNet.alpha;
        else
            warning('You did not specify a value for the relative penalization parameter (options.ElasticNet.alpha)')
            warning('I am using 0.5')
            ElasticNet_alpha = 0.5;
        end
    end
end

%********************************************************
%* Export every variable this function may have set.
%*
%* Two groups, deliberately written out longhand rather than via a
%* name-list + eval(): the unconditional group is always assigned above
%* (by a default if nothing else), so it's exported directly. The
%* conditional group is only assigned inside specific options.* branches,
%* so each is guarded by exist(...,'var') -- this is the one place we
%* genuinely don't know in advance what was set, and exist() answers that
%* without needing eval() to do it.
%********************************************************

% --- always set ---
opt.ny                       = ny;
opt.nunits                   = nunits;
opt.isMatlab                 = isMatlab;
opt.K                        = K;
opt.hor                      = hor;
opt.fhor                     = fhor;
opt.nethor                   = nethor;
opt.firstobs                 = firstobs;
opt.presample                = presample;
opt.noconstant                = noconstant;
opt.timetrend                = timetrend;
opt.minn_prior_tau           = minn_prior_tau;
opt.minn_prior_decay         = minn_prior_decay;
opt.minn_prior_lambda        = minn_prior_lambda;
opt.minn_prior_mu            = minn_prior_mu;
opt.minn_prior_omega         = minn_prior_omega;
opt.long_run_irf             = long_run_irf;
opt.irf_1STD                 = irf_1STD;
opt.cfrcst_yes               = cfrcst_yes;
opt.non_explosive_           = non_explosive_;
opt.heterosked                = heterosked;
opt.signs_irf                = signs_irf;
opt.narrative_signs_irf      = narrative_signs_irf;
opt.zeros_signs_irf          = zeros_signs_irf;
opt.proxy_irf                = proxy_irf;
opt.heterosked_irf           = heterosked_irf;
opt.hmoments_signs_irf       = hmoments_signs_irf;
opt.hmoments_eig_irf         = hmoments_eig_irf;
opt.noprint                  = noprint;
opt.nexogenous               = nexogenous;
opt.exogenous                = exogenous;
opt.cnnctdnss_               = cnnctdnss_;
opt.Ridge_                   = Ridge_;
opt.Lasso_                   = Lasso_;
opt.ElasticNet_               = ElasticNet_;
opt.set_irf                  = set_irf;
opt.robust_bayes_            = robust_bayes_;
opt.robust_credible_regions_ = robust_credible_regions_;
opt.exogenous_block          = exogenous_block;
opt.nz                       = nz;
opt.mixed_freq_on            = mixed_freq_on;
opt.dummy                    = dummy;
opt.flat                      = flat;
opt.priors                   = priors;
opt.varnames                 = varnames;
opt.y                        = y;   % possibly interpolated -- see header note

% --- only set under specific options.* branches ---
if exist('ww', 'var'),                  opt.ww                  = ww;                  end
if exist('K_shrinkage', 'var'),         opt.K_shrinkage          = K_shrinkage;         end
if exist('heterosked_regimes', 'var'),  opt.heterosked_regimes   = heterosked_regimes;  end
if exist('signs', 'var'),               opt.signs                = signs;               end
if exist('narrative', 'var'),           opt.narrative            = narrative;           end
if exist('zeros_signs', 'var'),         opt.zeros_signs          = zeros_signs;         end
if exist('f', 'var'),                   opt.f                    = f;                   end
if exist('sr', 'var'),                  opt.sr                   = sr;                  end
if exist('var_pos', 'var'),             opt.var_pos              = var_pos;             end
if exist('in', 'var'),                  opt.in                   = in;                  end
if exist('inols', 'var'),               opt.inols                = inols;               end
if exist('hmoments', 'var'),            opt.hmoments             = hmoments;            end
if exist('moment', 'var'),              opt.moment               = moment;              end
if exist('endo_path', 'var'),           opt.endo_path            = endo_path;           end
if exist('endo_path_index', 'var'),     opt.endo_path_index      = endo_path_index;     end
if exist('exo_index', 'var'),           opt.exo_index            = exo_index;           end
if exist('Omega', 'var'),               opt.Omega                = Omega;               end
if exist('L', 'var'),                   opt.L                    = L;                   end
if exist('opt_GiacomoniKitagawa', 'var'), opt.opt_GiacomoniKitagawa = opt_GiacomoniKitagawa; end
if exist('Ridge_lambda', 'var'),        opt.Ridge_lambda         = Ridge_lambda;        end
if exist('Lasso_lambda', 'var'),        opt.Lasso_lambda         = Lasso_lambda;        end
if exist('ElasticNet_lambda', 'var'),   opt.ElasticNet_lambda    = ElasticNet_lambda;   end
if exist('ElasticNet_alpha', 'var'),    opt.ElasticNet_alpha     = ElasticNet_alpha;    end
if exist('index_nan_var', 'var'),       opt.index_nan_var        = index_nan_var;       end
if exist('yoriginal', 'var'),           opt.yoriginal            = yoriginal;           end
if exist('yinterpol', 'var'),           opt.yinterpol            = yinterpol;           end
if exist('index', 'var'),               opt.index                = index;               end
if exist('prior', 'var'),               opt.prior                = prior;               end

end
