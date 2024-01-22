function [BVAR] = bvar_(y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bvar_' generates draws from the paramters of VAR model with error
% distributed as a multivariate normal

% Core Inputs:
% - y, data columns variables
% - lags, lag order of the VAR

% Additonal Inputs collected options:
% - options are not mandatory. if nothing is specified then a Jeffrey prior is
%   assumed
% - options.K, number of draws from the posterior distribution
% - options.hor, horizon to compute the impulse response
% - options.fhor, horizon to compute the out-of-sample forecast
% - options.priors, a string with he priors for the autoregressive paramters and
%   for the scaling matrix.
% (...) see below
% See the Hitchhiker's guide for more details. 
% https://github.com/naffe15/BVAR_/blob/master/HitchhikerGuide_.pdf


% Output: Draws from the conditional distribution of Phi, Sigma and
% Omega, impulse response with the various identification restrictions,
% forecast and marginal likelihood. 

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
% Revised, 27/02/2020
% Revised, 27/04/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('the bvar_ funtion needs at least two inputs: data and number of lags');
end
if lags < 1
    error('lags cannot be zero or negative');
end
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

% for mixed frequecy / irregurerly sampled data.
% Interpolate the missing values of each times series.
mixed_freq_on = 0;
if any(any(isnan(y))) %== 1
    if any(sum(isnan(y),2) == ny)
        warning('Cannot activate the Mixed Frequency BVAR. Data contains raws that have only ''nan''s. Each variable is interpolated individually.')        
    else
        disp('Activating the Mixed Frequency BVAR')
        mixed_freq_on = 1;
    end    
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
priors   = priors_( );

% declaring the names for the observable variables
for v = 1 : ny
    eval(['varnames{'   num2str(v) '} =  ''Var' num2str(v) ''';'])
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
        [f,sr] = sign2matrix(zeros_signs,ny);
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
%* Consistency Checks
%********************************************************
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
if dummy == 1 && nexogenous > 0 
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
    %fprintf('To change this, include the exogenous forecasts in options.exogenous.\n')
    exogenous = [exogenous; zeros(fhor,nexogenous)];
end

%********************************************************
%* Priors and Posterior Distributions
%********************************************************

idx = firstobs+presample-lags:firstobs+nobs-1;
nx  = 1;
if noconstant
    nx = 0;
end

% organize data as  yy = XX B + E
if nunits == 1
    [yy,XX] = YXB_(y(idx, :),lags,[nx timetrend]);
else % pooled units
    yy = []; XX = [];
    for nunit = 1 : nunits
        [yy,XX] = YXB_(y(idx, :, nunit),lags,[nx timetrend]);
        XX = [XX; XX];
        yy = [yy; yy];
    end
end

ydum    = [];
xdum    = [];
pbreaks = 0;
lambda  = 0;
mu      = 0;

ydata   = y(idx, :, :);
T       = size(ydata, 1);
if T-lags < lags*ny + nx %+ flat*(ny+1)
    error('Less observations than regressors: increase the # of obs or decrease the # of lags.')
end
xdata   = ones(T,nx);
if timetrend == 1
    % xdata = [xdata [1:T]'];
    xdata = [xdata [1-lags : T-lags]'];
end
if nexogenous > 0
    xdata = [xdata exogenous(idx,:)]; 
    XX    = [XX exogenous(idx(1)+lags : idx(end),:)];
end

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
    varols = fast_ols(varols.y,varols.X);
end

if Ridge_ == 1
    varRidge   = varols;
    varRidge.B = (eye(size(varRidge.B,1)) + Ridge_lambda*varols.xxi)\varRidge.B;
    % varRidge.B = inv(eye(size(varRidge.B,1)) + Ridge_lambda*varols.xxi)*varRidge.B;
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


% specify the prior
if dummy == 1
    % MINNESOTA PRIOR:
    mnprior.tight         = minn_prior_tau;
    mnprior.decay         = minn_prior_decay;
    % Use only initializations lags for the variance prior
    % vprior.sig            = std(y(firstobs+presample-lags : firstobs+presample,:))';
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
%     varp            = rfvar3([ydata(1:Tp, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
    varp           = rfvar3([y(firstobs-lags : firstobs+presample-1, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
    Tup            = size(varp.u, 1);
    prior.df       = Tup - ny*lags - nx - flat*(ny+1);
    prior.S        = varp.u' * varp.u;
    prior.XXi      = varp.xxi;
    prior.PhiHat   = varp.B;
    prior.YYdum    = varp.y;
    prior.XXdum    = varp.X;
    if prior.df < ny
        error('Too few degrees of freedom in the Inverse-Wishart part of prior distribution. You should increase training sample size.')
    end
    prior.minn_prior_tau    = minn_prior_tau;
    prior.minn_prior_decay  = minn_prior_decay;
    prior.minn_prior_lambda = minn_prior_lambda;
    prior.minn_prior_mu     = minn_prior_mu;
    prior.minn_prior_omega  = minn_prior_omega;
    
elseif dummy == 0
    % JEFFREY OR UNIFORMATIVE PRIOR:
    % [priors] = jeffrey(y,lags);
    prior.name  = 'Jeffrey';
    
elseif dummy == 2
    % MN-IW
    prior.name  = 'MultivariateNormal-InverseWishart';    
end

% specify the posterior (the varols agin on actual+dummy)
[posterior,var] = posterior_(y);

% compute the second and fourth moments robust to error distribution misspecification  
if robust_bayes_ > 0 
    
    
    %% Kurtosis
    %var.u
    Dpl       = duplicationmatrix(ny);
    Dplus     = pinv(Dpl);
    I_ny      = vech(eye(ny));
    %I__ny     = reshape(eye(ny),ny^2,1); 
    Sig_      = 1/nobs * varols.u' * varols.u;
    vech_Sig_ = vech(Sig_);

    % whithen the reduced form errors
    SigChol   = chol(Sig_,'lower');
    iSigChol  = inv(SigChol);
    white_u   = varols.u * iSigChol';     
    % compute fourth centered moments
    Khat_     = fourthmom(white_u);
    % shrink towards the normal multivariate fourth mom
    % equation (17) at page 161
    K1_       = nobs/(nobs+K_shrinkage) * (Khat_ - I_ny * I_ny') ; 
    K2_       = K_shrinkage/(nobs+K_shrinkage) * ...
        Dplus *(eye(ny^2) + commutationmatrix(ny)) * Dplus'; 
        %Dplus *(eye(ny^2) + commutationmatrix(ny) + I__ny*I__ny') * Dplus'; 
    Kstar_     = K1_ + K2_;
    % equation (8) of remark 3 at page 158
    Left  = Dplus * kron(SigChol,SigChol) * Dpl;
    Rigth = Left';
    
    vech_Sig_cov_            = 1/nobs * Left * (Kstar_) * Rigth;   
    vech_Sig_cov_lower_chol  = chol(vech_Sig_cov_,'lower');
    
    if robust_bayes_ > 1 
        %% Skewness
        ivech_Sig_cov_ = pinv(vech_Sig_cov_);
        Sstar_         = nobs /(nobs+K_shrinkage) * thirdmom(varols.u);
        mu_cov_        = 1/nobs * Sstar_ * ivech_Sig_cov_ * Sstar_';
        mu_cov_chol    = chol(Sig_ - mu_cov_,'lower');
    end
end

%*********************************************************
%* Compute the log marginal data density for the VAR model
%*********************************************************

try
    posterior_int = matrictint(posterior.S, posterior.df, posterior.XXi);
catch
    warning('I could not compute the marginal likelihood');
    posterior_int = nan;
end
        
if dummy == 1 % only for minnesota dummy
    prior_int = matrictint(prior.S, prior.df, prior.XXi);
    lik_nobs  = posterior.df - prior.df;
    log_dnsty = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);

elseif dummy == 0 % jeffrey prior
    lik_nobs  = posterior.df;
    log_dnsty = posterior_int - 0.5*ny*lik_nobs*log(2*pi);

elseif dummy == 2 % conjugate MN-IW prior
    log_dnsty = mniw_log_dnsty(prior,posterior,varols);
    
end


%**************************************************
%* Generating draws form the Posterior Distribution
%**************************************************

% Preallocation of memory
% Matrices for collecting draws from Posterior Density
% Last dimension corresponds to a specific draw
Phi_draws     = zeros(ny*lags+nx+timetrend + nexogenous, ny, K);   % Autoregressive Parameters
Sigma_draws   = zeros(ny,ny,K);                     % Shocks Covariance
ir_draws      = zeros(ny,hor,ny,K);                 % variable, horizon, shock and draws - Cholesky IRF
irlr_draws    = zeros(ny,hor,ny,K);                 % variable, horizon, shock and draws - Long Run IRF
Qlr_draws     = zeros(ny,ny,K);                     % long run impact matrix
e_draws       = zeros(size(yy,1), ny,K);                  % residuals
yhatfut_no_shocks         = NaN(fhor, ny, K, nunits);   % forecasts with shocks
yhatfut_with_shocks       = NaN(fhor, ny, K, nunits);   % forecast without the shocks
yhatfut_cfrcst            = NaN(fhor, ny, K, nunits);   % forecast conditional on endogenous path
logL                      = NaN(K,1);
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
    % Ox_draws  = Sigma_draws;
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
forecast_data.initval     = ydata(end-lags+1:end, :, :);

% Settings for the MFVAR
if mixed_freq_on == 1
    KFoptions.index   = index;
    KFoptions.noprint = noprint;
end


try
    S_inv_upper_chol    = chol(inv(posterior.S));
catch
    warning('POSTERIOR MEAN of SIGMA IS ILL-BEHAVED (NON-POSITIVE DEFINITE)')
    try
        warning('I will try with the LDL decomposition')
        iS               = inv(posterior.S);
        [Left, D, Right]        = ldl(iS);
        S_inv_upper_chol = sqrt(D) * Left' * Right';%chol()
    catch
        warning('I will add +1e-05 to the diagonal')
        S_inv_upper_chol    = chol(inv(posterior.S + 1e-05*eye(ny)));
    end
end


XXi_lower_chol      = chol(posterior.XXi)';

nk                  = ny*lags+nx+timetrend + nexogenous;

% Declaration of the companion matrix
Companion_matrix = diag(ones(ny*(lags-1),1),-ny);
p  = 0;
dd = 0;

waitbar_yes = 0;
if K > 99
    waitbar_yes = 1;
    wb = waitbar(0, 'Generating draws from the Posterior Distribution');
end

for  d =  1 : K
    
    while dd == 0 
        %======================================================================
        % Inferece: Drawing from the posterior distribution
        % Step 1: draw from the Covariance
        if robust_bayes_ > 0 % robust to kurtosis 
            tec = 0;
            while tec == 0
                Sig0  = randn((ny+1)*ny/2, 1);
                Sig1  = vech_Sig_cov_lower_chol * Sig0;
                Sig2  = vech_Sig_ + Sig1;
                Sigma = ivech(Sig2);
                %vech_Sig_ = vech(Sig_);
                %vech_Sig_cov_lower_chol  = 1/nobs * (Kstar_ - vech_Sig_ * vech_Sig_');
                [Sigma_lower_chol,flag] = chol(Sigma,'lower');
                if flag == 0
                    tec = 1;
                end
            end
        else    
            Sigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);
            Sigma_lower_chol = chol(Sigma)';
        end

        Phi1 = randn(nk * ny, 1);
        Phi2 = kron(Sigma_lower_chol , XXi_lower_chol) * Phi1;
        Phi3 = reshape(Phi2, nk, ny);
        Phi  = Phi3 + posterior.PhiHat;
        
        if robust_bayes_ > 1 % robust to kurtosis and skewness
            mu_rob           = posterior.PhiHat(ny*lags+1,:) + (Sstar_* ivech_Sig_cov_ * (Sig2  - vech_Sig_))';
            % as if assuming no lags (only costant) (X'X) = T
            % mu_cov_chol      = chol(Sigma - mu_cov_,'lower');            
            Phi(ny*lags+1,:) = mu_rob + (mu_cov_chol*randn(ny, 1))';
        end
        
        Companion_matrix(1:ny,:) = Phi(1:ny*lags,:)';
        test = (abs(eig(Companion_matrix)));
        if non_explosive_ == 1
            % Checking the eigenvalues of the companion matrix (on or inside the
            % unit circle)
            if any(test>1.01) %any(test>1.0000000000001)
                p = p+1;                
            else
                dd = 1;
            end
        else
            dd = 1;
        end
    end
    
    % store the draws
    Phi_draws(:,:,d)   = Phi;
    Sigma_draws(:,:,d) = Sigma;
    errors             = yy - XX * Phi;
    e_draws(:,:,d)     = errors;
    
    
    %======================================================================
    % IRF
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
                    %waitbar(ell/L, wb1);
                end
                % irsign_upper_draws(:,:,:,d) = max(irsign0,[],4,''omitnan'');
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
                    [irsign0(:,:,:,ell), ~] = iresponse_sign(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,narrative);
                    %waitbar(ell/L, wb1);
                end
                % irsign_upper_draws(:,:,:,d) = max(irsign0,[],4,''omitnan'');
                irnarrsign_upper_draws(:,:,:,d) = max(irsign0,[],4,"omitnan");
                irnarrsign_lower_draws(:,:,:,d) = min(irsign0,[],4,"omitnan");
            end
        end

    end
    % with zeros and sign restrictions
    if zeros_signs_irf == 1         %= iresponse_zeros_signs( Phi,Sigma,bvar1.hor,lags,var_pos,f,sr);
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
                    %waitbar(ell/L, wb1);
                end
                % irsign_upper_draws(:,:,:,d) = max(irsign0,[],4,''omitnan'');
                irzerosign_upper_draws(:,:,:,d) = max(irsign0,[],4,"omitnan");
                irzerosign_lower_draws(:,:,:,d) = min(irsign0,[],4,"omitnan");
            end
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
        Omegap_draws(:,:,d) = icholSig*[[omega11; omega21'] , [omega12; omega22]];
        %max(max(abs(Omegap_draws(:,:,d)*Omegap_draws(:,:,d)'-Sigma)))
        % max(max(abs(Omegap_draws(:,:,d)*Omegap_draws(:,:,d)'-eye(ny))))
        clear tmp_
    end
    % with heteroskedasticity 
    if heterosked_irf == 1
        [irheterosked,Omegah]       = iresponse_heterosked(Phi(1 : ny*lags, 1 : ny),errors,hor,heterosked_regimes);
        irheterosked_draws(:,:,:,d) = irheterosked;
        Omegah_draws(:,:,d)         = Omegah;
    end
    % with higher-moments and sign restrictions
    if hmoments_signs_irf == 1
        [irhmomsign,Omega]         = iresponse_sign_hmoments(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,hmoments,f);
        irhmomsign_draws(:,:,:,d)  = irhmomsign;
        Omegam_draws(:,:,d)        = Omega;
%         if isnan(Omega)
%             nan_count = nan_count + 1;
%             disp([nan_count nan_count/d])
%         end
    end
    % with higher-moments eigenvalue decomposition
    if hmoments_eig_irf == 1
        [irhmomeig,Omega]         = iresponse_sign_hmoments_eig(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,moment);
        irhmomeig_draws(:,:,:,d)  = irhmomeig;
        Omegae_draws(:,:,d)        = Omega;
    end
    
    %======================================================================
    % Forecasts
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
    % Mixed Frequency
    if mixed_freq_on
        % Checking the eigenvalues of the companion matrix (on or inside the
        % unit circle). Needed for the initial of KF
        Companion_matrix(1:ny,:) = Phi(1:ny*lags,:)';
        %test = (abs(eig(Companion_matrix)));
        if any(test>1.0000000000001)
            KFoptions.initialCond = 1;
        end
        % Forward Kalman Filter and Smoother
        [logL,KFout]     = kfilternan(Phi,Sigma,yoriginal,KFoptions);
        yfill(:,:,d)     = KFout.smoothSt_plus_ss(:,KFout.index_var);
        yfilt(:,:,d)     = KFout.filteredSt_plus_ss(:,KFout.index_var);
        % recompute the posterior with smoothed data
        [posterior1]     = posterior_(yfill(:,:,d));
        S_inv_upper_chol = chol(inv(posterior1.S));
        XXi_lower_chol   = chol(posterior1.XXi)';
        posterior.PhiHat = posterior1.PhiHat;
        forecast_data.initval = yfill(end-lags+1:end, :, d);
        logL(d) = KFout.logL;
    end
    
    %======================================================================
    % Connectedness
    if cnnctdnss_ > 0 
%         if use_omega == 1
%             [C] = connectedness(Phi,Sigma,nethor,Omega); 
%         else
%         end
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
    
    
    if waitbar_yes, waitbar(d/K, wb); end
    dd = 0; % reset 
end
if waitbar_yes, close(wb); end


%********************************************************
%* Storing the resutls
%*******************************************************

%==========================================================================
% classical inference: OLS estimator
BVAR.Phi_ols    = varols.B;
BVAR.e_ols      = varols.u;
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
        %[H,Pv] = kstest(BVAR.e_ols(:,gg)/sqrt(BVAR.Sigma_ols(gg,gg)));
        [H,Pv] = kstest(BVAR.white_e_ols(:,gg));
        BVAR.HP(gg,:) = [H,Pv];
        %      H = 0 => Do not reject the null hypothesis at the 5% significance
        %      level. 
    end
else
    BVAR.HP = [];
end
%==========================================================================
% Penalized Approaches (Regularization)
% Parameters Estiamtes
% if Ridge_ == 1
%     BVAR.Ridge.Phi    = varRidge.B;
%     BVAR.Ridge.e      = varRidge.u;
%     BVAR.Ridge.Sigma  = 1/(nobs-nk)*varRidge.u'*varRidge.u;
%     % info crit
%     [BVAR.Ridge.InfoCrit.AIC, BVAR.Ridge.InfoCrit.HQIC, BVAR.Ridge.InfoCrit.BIC] ...
%         = IC(BVAR.Ridge.Sigma, BVAR.Ridge.e, nobs, nk);
%     % ir with recursive identification
%     BVAR.Ridge.ir      = iresponse(BVAR.Ridge.Phi,BVAR.Ridge.Sigma,hor,eye(ny));
%     % connectedness 
% end
% if Lasso_ == 1
%     BVAR.Lasso.Phi    = varLasso.B;
%     BVAR.Lasso.e      = varLasso.u;
%     BVAR.Lasso.Sigma  = 1/(nobs-nk) * varLasso.u'*varLasso.u;
%     % info crit
%     [BVAR.Lasso.InfoCrit.AIC, BVAR.Lasso.InfoCrit.HQIC, BVAR.Lasso.InfoCrit.BIC] ...
%         = IC(BVAR.Lasso.Sigma, BVAR.Lasso.e, nobs, nk);
%     % ir with recursive identification
%     BVAR.Lasso.ir      = iresponse(BVAR.Lasso.Phi,BVAR.Lasso.Sigma,hor,eye(ny));   
% end
% if ElasticNet_ == 1
%     BVAR.ElasticNet.Phi    = varElasticNet.B;
%     BVAR.ElasticNet.e      = varElasticNet.u;
%     BVAR.ElasticNet.Sigma  = 1/(nobs-nk) * varElasticNet.u'*varElasticNet.u;
%     % info crit
%     [BVAR.ElasticNet.InfoCrit.AIC, BVAR.ElasticNet.InfoCrit.HQIC, BVAR.ElasticNet.InfoCrit.BIC] ...
%         = IC(BVAR.ElasticNet.Sigma, BVAR.ElasticNet.e, nobs, nk);
%     % ir with recursive identification
%     BVAR.ElasticNet.ir      = iresponse(BVAR.ElasticNet.Phi,BVAR.ElasticNet.Sigma,hor,eye(ny));   
% end
% Parameters Estiamtes - IRFs - Connectedness
%set_irf = signs_irf + narrative_signs_irf + zeros_signs_irf;
penalizationstrn = {'Ridge','Lasso','ElasticNet'};
% loop across penalization approaches
for pp = 1 : 3    
    Phi_ = []; Sigma_ = []; u_ = []; Omega_(:,:,1) = eye(ny);  InfoCrit_ = [];   
    if eval([ penalizationstrn{pp} '_ == 1'])
        % AR param
        eval(['BVAR.' penalizationstrn{pp} '.Phi    = var' penalizationstrn{pp} '.B;']);
        eval(['Phi_   = BVAR.' penalizationstrn{pp} '.Phi;'])       
        % Error term
        eval(['BVAR.' penalizationstrn{pp} '.e      = var' penalizationstrn{pp} '.u;']);
        eval(['u_ = BVAR.' penalizationstrn{pp} '.e;'])
        % Covariance Matrix        
        eval(['BVAR.' penalizationstrn{pp} '.Sigma  = 1/(nobs-nk) * var' penalizationstrn{pp} '.u'' * var' penalizationstrn{pp} '.u;'])        
        eval(['Sigma_ = BVAR.' penalizationstrn{pp} '.Sigma;'])
        % Recursive IRFs
        eval(['BVAR.' penalizationstrn{pp} '.ir      = iresponse(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,eye(ny));']);                
        % info crit
        [InfoCrit_.AIC, InfoCrit_.HQIC, InfoCrit_.BIC] ...
            = IC(Sigma_, u_, nobs, nk);
        eval(['BVAR.' penalizationstrn{pp} '.InfoCrit      = InfoCrit_;']);                                     
        % IRFs with different identification schemes
        % point identification: LR
        if long_run_irf == 1
            eval(['[BVAR.' penalizationstrn{pp} '.irlr,Omega_] = iresponse_longrun(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,lags);'])
        end
        % point identification: with heteroskedasticity
        if heterosked_irf == 1
            eval(['[BVAR.' penalizationstrn{pp} '.irheterosked,Omega_] = iresponse_heterosked(Phi_(1 : ny*lags, 1 : ny),u_,hor,heterosked_regimes);'])
        end
        % set identification: signs
        if set_irf > 0  
            wb = waitbar(0, ['Generating rotations for set-identification - ' penalizationstrn{pp} ' Estimator']);
        end        
        if signs_irf == 1
            for d1 = 1 : set_irf                
                eval(['[BVAR.' penalizationstrn{pp} '.irsign_draws(:,:,:,' num2str(d1) ...
                    '),Omega_(:,:,' num2str(d1) ')] = iresponse_sign(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,signs);'])
                waitbar(d1/set_irf, wb);
            end            
        end
        % set identification: narrative and sign restrictions
        if narrative_signs_irf == 1
            for d1 = 1 : set_irf
                eval(['[BVAR.' penalizationstrn{pp} '.irnarrsign_draws(:,:,:,' num2str(d1) ...
                    '),Omega_(:,:,' num2str(d1) ')] = iresponse_sign_narrative(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,signs,narrative);'])
                waitbar(d1/set_irf, wb);
            end
        end
        % set identification: zeros and sign restrictions
        if zeros_signs_irf == 1
            for d1 = 1 : set_irf
                eval(['[BVAR.' penalizationstrn{pp} 'irzerosign_draws.(:,:,:,' num2str(d1) ...
                    '),Omega_(:,:,' num2str(d1) ')] = iresponse_zeros_signs(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,lags,var_pos,f,sr);'])
                waitbar(d1/set_irf, wb);
            end
        end
        if set_irf>0, close(wb); end
    end
    if cnnctdnss_ == 1 % default identification (Pesaran and Shin)
        eval(['BVAR.' penalizationstrn{pp} '.Connectedness = connectedness(Phi_(1 : ny*lags, 1 : ny),Sigma_,nethor);']);        
    elseif cnnctdnss_ == 2 % customized identification
        Sigma_lower_chol = chol(Sigma_)';
        Omegam           = median(Omega_,3);
        eval(['BVAR.' penalizationstrn{pp} '.Connectedness = connectedness(Phi_(1 : ny*lags, 1 : ny), Sigma_, nethor, Sigma_lower_chol * Omegam);']);                
    end

end
    
%==========================================================================
% bayesian inference:
% the last dimension of these objects corresponds to a draw from the posterior

% inference and IRFs
BVAR.Phi_draws    = Phi_draws;          % draws from the autoregressive part
BVAR.Sigma_draws  = Sigma_draws;        % draws from the covarance matrix
BVAR.alpha_draws  = Phi_draws;          % Older name: draws from the autoregressive part
BVAR.sigma_draws  = Sigma_draws;        % Older name: draws from the covarance matrix

BVAR.ir_draws     = ir_draws;           % draws from the IRF with cholesky
BVAR.irlr_draws   = irlr_draws;         % draws from the IRF with Long Run
BVAR.Qlr_draws    = Qlr_draws;          % Long Run Rotation matrix
BVAR.lags         = lags;               % lags
BVAR.N            = ny;                 % number of variables
BVAR.e_draws      = e_draws;            % residuals
BVAR.e            = e_draws;            % backward compatible with earlier versions

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
end
%
BVAR.varnames     = varnames;
BVAR.ndraws       = K;

if signs_irf == 1 && narrative_signs_irf == 0
    BVAR.irsign_draws = irsign_draws;
    BVAR.Omegas       = Omega_draws;
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    %fprintf('\nPosterior plausibility of identifying restrictions: %0.4g.\n',BVAR.postPlaus)
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irsign_lower_draws = irsign_lower_draws;
        BVAR.irsign_upper_draws = irsign_upper_draws;
        for ss = 1 : ny
            % flip the order
            % from: variable, horizon, [shock], draw 
            % to:   draw    , horizon, [shock], variable
            rMinPost = permute(squeeze(irsign_lower_draws(:,:,ss,:)), [3 2 1] );
            rMaxPost = permute(squeeze(irsign_upper_draws(:,:,ss,:)), [3 2 1] );
            % Compute robustified credible region for IS. 
            % [integrate out draws]
            [credlb,credub] = credibleRegion(rMinPost,rMaxPost,opt_GiacomoniKitagawa);
            BVAR.irsign_robust_credible_bands.u(:,:,ss) = credub';
            BVAR.irsign_robust_credible_bands.l(:,:,ss) = credlb';

            % % Compute highest posterior density (HPD) interval under single prior.
            % [hpdlb,hpdub] = highestPosteriorDensity(rSinglePriorPost,opt);
            % 
            % postMeanBoundWidth = meanub - meanlb; % Width of posterior mean bounds.
            % hpdWidth = hpdub - hpdlb; % Width of highest posterior density regions
            % credWidth = credub - credlb; % Width of robustified credible region
            % 
            % priorInformativeness = 1 - hpdWidth./credWidth; % Informativeness of prior
        end
    end
else
    BVAR.irsign_draws = [];
    BVAR.Omegas       = [];
end
if narrative_signs_irf == 1
    BVAR.irnarrsign_draws = irnarrsign_draws;
    BVAR.Omegan           = Omegan_draws;
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    %fprintf('\nPosterior plausibility of identifying restrictions: %0.4g.\n',BVAR.postPlaus)
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irnarrsign_lower_draws = irnarrsign_lower_draws;
        BVAR.irnarrsign_upper_draws = irnarrsign_upper_draws;
        for ss = 1 : ny
            % flip the order
            % from: variable, horizon, [shock], draw 
            % to:   draw    , horizon, [shock], variable
            rMinPost = permute(squeeze(irzerosign_lower_draws(:,:,ss,:)), [3 2 1] );
            rMaxPost = permute(squeeze(irzerosign_upper_draws(:,:,ss,:)), [3 2 1] );
            % Compute robustified credible region for IS. 
            % [integrate out draws]
            [credlb,credub] = credibleRegion(rMinPost,rMaxPost,opt_GiacomoniKitagawa);
            BVAR.irnarrsign_robust_credible_bands.u(:,:,ss) = credub';
            BVAR.irnarrsign_robust_credible_bands.l(:,:,ss) = credlb';
        end
    end    
else
    BVAR.irnarrsign_draws = [];
    BVAR.Omegan           = [];
end
if zeros_signs_irf == 1
    BVAR.irzerosign_draws   = irzerosign_draws;
    BVAR.Omegaz             = Omegaz_draws;
    BVAR.postPlaus          = (1 - sum(OmegaEmpty)/K);
    %fprintf('\nPosterior plausibility of identifying restrictions: %0.4g.\n',BVAR.postPlaus)
    if robust_credible_regions_
        BVAR.robust_credible_regions_ = robust_credible_regions_;
        BVAR.irzerosign_lower_draws = irzerosign_lower_draws;
        BVAR.irzerosign_upper_draws = irzerosign_upper_draws;
        for ss = 1 : ny
            % flip the order
            % from: variable, horizon, [shock], draw 
            % to:   draw    , horizon, [shock], variable
            rMinPost = permute(squeeze(irzerosign_lower_draws(:,:,ss,:)), [3 2 1] );
            rMaxPost = permute(squeeze(irzerosign_upper_draws(:,:,ss,:)), [3 2 1] );
            % Compute robustified credible region for IS. 
            % [integrate out draws]
            [credlb,credub] = credibleRegion(rMinPost,rMaxPost,opt_GiacomoniKitagawa);
            BVAR.irzerosign_robust_credible_bands.u(:,:,ss) = credub';
            BVAR.irzerosign_robust_credible_bands.l(:,:,ss) = credlb';
        end
    end
else
    BVAR.irzerosign_draws = [];
    BVAR.Omegaz           = [];
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
end

% bvar connetedness
if cnnctdnss_
    BVAR.Connectedness.Index         = CnndtnssIndex;
    BVAR.Connectedness.FromAlltoUnit = CnndtnssFromAlltoUnit;
    BVAR.Connectedness.FromUnitToAll = CnndtnssFromUnitToAll;
    BVAR.Connectedness.theta         = Ctheta;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of bvar_.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%********************************************************
%********************************************************
    function [posterior,var] = posterior_(y)
        
        % This part is needed in the case of missing values. The posterior
        % needs to be revaluated given a new value of the Kalman smoothed
        % observables
        %===========================================
        ydata = y(idx, :, :);
        %===========================================
        T     =  size(ydata, 1);
        xdata = ones(T,nx);
        if timetrend ==1
            xdata = [xdata [1-lags:T-lags]'];
            % xdata = [xdata [1:T]'];
        end
        if nexogenous >0
            xdata = [xdata exogenous(idx,:)];
        end
        % posterior density
        if nunits == 1
            if heterosked == 0
                var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu);
            else
                var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu, ww);
            end
        else  % pooled units
            var.y = []; var.X = [];
            for nunt = 1 : nunits
                tmp_var = rfvar3(ydata(:,:,nunt), lags, xdata, [T; T], 0, 0);
                var.X = [var.X; tmp_var.X];
                var.y = [var.y; tmp_var.y];
            end
            % Compute OLS regression and residuals of the pooled estimator
            var = fast_ols(var.y,var.X);
        end
        Tu = size(var.u, 1);
                
        if dummy ==  1
            %********************************************************
            % Minnesota Prior
            %********************************************************
            % Prior density
            Tp = presample + lags;
            if nx
                xdata = xdata(1:Tp, :);
            else
                xdata = [];
            end
            % varp            = rfvar3([ydata(1:Tp, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
            varp            = rfvar3([y(firstobs-lags : firstobs+presample-1, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
            Tup             = size(varp.u, 1);
            prior.df        = Tup - ny*lags - nx - flat*(ny+1);
            prior.S         = varp.u' * varp.u;
            prior.XXi       = varp.xxi;
            prior.PhiHat    = varp.B;
            priors.YYdum    = varp.y;
            priors.XXdum    = varp.X;
            if prior.df < ny
                error('Too few degrees of freedom in the Inverse-Wishart part of prior distribution. You should increase training sample size.')
            end
            posterior.df    = Tu - ny*lags - nx - flat*(ny+1);
            posterior.S     = var.u' * var.u;
            posterior.XXi   = var.xxi;
            posterior.PhiHat = var.B;
            
        elseif dummy == 2
            %********************************************************
            % Conjugate Prior
            %********************************************************
            Ai              = inv(prior.Phi.cov);
            posterior.df    = Tu - ny*lags - nx - nexogenous - prior.Sigma.df;
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
        end
        
    end
%********************************************************
%********************************************************
    function [priors] = priors_( )
        
        priors.name  = 'N/A';

    end
%********************************************************
%********************************************************
    function out = fast_ols(y,X)
        % Compute OLS regression and residuals
        [vl,d_,vr] = svd(X,0);
        di = 1./diag(d_);
        B = (vr.*repmat(di',ny*lags+nx,1))*vl'*y;
        u = y-X*B;
        xxi = vr.*repmat(di',ny*lags+nx,1);
        xxi = xxi*xxi';

        out.B = B;
        out.u = u;
        out.xxi = xxi;
        out.y   = y;
        out.X   = X;
    end
%********************************************************
%********************************************************

end
