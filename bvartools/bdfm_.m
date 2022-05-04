function [BDFM] = bdfm_(y,lags,nfac,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bdfm_' generates draws from the paramters and latent variables of
% static/dynamic factor model with factors and idiosyncrati errors
% distributed as multivariate normals 

% Core Inputs:
% - y, data columns variables
% - lags, lag order of the dynamic factor model (when =0 static factor model)
% - nfac, number of factors 

% Additonal Inputs collected options:
% - options are not mandatory. if nothing is specified then a number of
% default options are assumed. 
% (...) see below
% See the Hitchhiker's guide for more details. 
% https://github.com/naffe15/BVAR_/blob/master/HitchhikerGuide_.pdf

% Output: Draws from the conditional distribution of Phi, Sigma and
% Omega, impulse response with the cholesky decompotion and long run
% restricions, forecast and marginal likelihood.

% Filippo Ferroni, 7/1/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    error('the BDFM function needs at least three inputs: data, # of lags and # of factors');
end
if lags < 0
    error('lags cannot be negative');
end
if nfac < 1 
    error('need at least one factor');
end

% length of TS and number of observable variables
[T,ny] = size(y);

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

%********************************************************
%* DEFAULT SETTINGS
%********************************************************
%
ilags               = 0;            % lag strucutre of the idiosyncratic component
K                   = 1000;        % number of draws from the posterior
hor                 = 24;           % horizon for the IRF
fhor                = 12;           % horizon for the forecasts
noconstant          = 1;            % when 0, includes a constatn in the factor VAR
l_F                 = 0;            % factors entering with lags in the measurement equation
l_G                 = 0;            % measurement errors entering with lags in the measurement equation
sig_fix             = [];           % when non empty, fix the variance of the IID part in the idiosyncratic error (if empty not-fixed)
noise_var           = 0.2;          % size of the noise variance when perturbing the PC factors to initiallize the Gibbs sampler
nexogenous          = 0;            % exogenous variables in the factor equations
timetrend           = 0;            % time trend in the factor equations
long_run_irf        = 0;            % when 0, it does not compute long run IRF
signs_irf           = 0;
narrative_signs_irf = 0;
proxy_irf           = 0;
do_the_demean       = 1;            % demean the data
dummy               = 2;            % Conjugate MN-IW on the factors
use_the_ss_kalman_gain = 0;
burnin_             = 5000;
skip_               = 20;
% Default Prior Distributions 
% Priors for parameters related to factors (F)
% Lambda
priors.F.Lambda.mean = 0;
priors.F.Lambda.cov = 10;
% Phi
priors.F.Phi.mean                  = zeros(nfac*lags + (1-noconstant) + timetrend + nexogenous ,  nfac);
priors.F.Phi.mean(1:nfac , 1:nfac) = 0.2 * ones(nfac);
priors.F.Phi.cov                   = 5   * eye((nfac*lags+(1-noconstant) + timetrend + nexogenous));
% Sigma
priors.F.Sigma.scale = eye(nfac); 
priors.F.Sigma.df    = nfac + nexogenous + timetrend + 1;

% Priors for parameters related to idiosyncratic term (G)
% phi
priors.G.Phi.mean = 0.0;
priors.G.Phi.cov  = 5.0;
% sigma
priors.G.Sigma.scale = 1;
priors.G.Sigma.df    = 4;


% declaring the names for the observable variables
for v = 1 : ny
    eval(['varnames{'   num2str(v) '} =  ''Var' num2str(v) ''';'])
end

%********************************************************
%* CUSTOMIZED SETTINGS
%********************************************************
if nargin > 3
    if isfield(options,'vnames')==1
        varnames = options.vnames;
    end
    %======================================================================
    % Inference options
    %======================================================================
    if isfield(options,'K')==1
        % number of draws from the posterior distribution
        K = options.K;
    end
    if isfield(options,'burnin')==1
        % number of burnin draws from the posterior distribution
        burnin_ = options.burnin;        
    end
    if isfield(options,'skip')==1
        % keep one ecvery skip_ draws
        skip_ = options.skip;        
    end
    if isfield(options,'noconstant')==1
        % when 0, includes a constatn in the factor VAR
        noconstant = options.noconstant;
    end
    if isfield(options,'ilags')==1
        % lag structure of the idiosyncratic component
        ilags 	  =	options.ilags;
        if ilags > 0
            warning(['AR1 for the idiosyncratic errors is currently not supported by the toolbox. '...
                'Assume that all AR1 coefficients of the idiosyncratic errors are zero.'])
            ilags = 0;
        end
    end    
    if isfield(options,'l_F')==1
        % factors entering with lags in the measurement equation
        l_F       =	options.l_F;
    end
    if isfield(options,'l_G')==1
        % measurement errors entering with lags in the measurement equation
        l_G       =	options.l_G;
    end
    if isfield(options,'sig_fix')==1
        % when sig_fix is not-empty, fix the variance of the IID part in the idiosyncratic error (if empty not-fixed)
        sig_fix   =	options.sig_fix	;
    end
    if isfield(options,'noise_var')==1        
        % size of the noise variance when perturbing the PC factors to initiallize the Gibbs sampler
        noise_var =	options.noise_var;
    end
    if isfield(options,'do_the_demean')==1
        % demean the data unless otherwise specified (do_the_demean=0)
        do_the_demean = options.do_the_demean;
    end
    if isfield(options,'use_the_ss_kalman_gain')==1
        use_the_ss_kalman_gain = options.use_the_ss_kalman_gain;
    end
    %======================================================================
    % Conjugate/Hierachical MN-IW prior options
    %======================================================================
%     if (isfield(options,'priors')==1 && strcmp(options.priors.name,'Conjugate')==1) || (isfield(options,'priors')==1 && strcmp(options.priors.name,'conjugate')==1) || ...
%        (isfield(options,'prior')==1 && strcmp(options.prior.name,'Conjugate')==1) || (isfield(options,'prior')==1 && strcmp(options.prior.name,'conjugate')==1)
    if isfield(options,'priors')==1  || isfield(options,'priors')==1 
        
        if isfield(options,'prior')==1
            options.priors = options.prior;
        end
        priors.name= 'Conjugate';
        
        % =========================================================== %
        % Priors for the factors
        % =========================================================== %
        if isfield(options.priors,'F') == 1  
            
            % =========================================================== %
            % Priors for the Factor AR parameters
            % Default values
            %               priors.F.Phi.mean                  = zeros(nfac*lags + (1-noconstant) + timetrend + nexogenous ,  nfac);
            %               priors.F.Phi.mean(1:nfac , 1:nfac) = 0.2 * ones(nfac);
            %               priors.F.Phi.cov                   = 5   * eye((nfac*lags+(1-noconstant) + timetrend + nexogenous));
            if isfield(options.priors.F,'Phi') == 1
                % mean
                if isfield(options.priors.F.Phi,'mean') == 1
                    priors.F.Phi.mean  = options.priors.F.Phi.mean;
                    if max(size(priors.F.Phi.mean) ~= [nfac*lags+(1-noconstant)+timetrend+nexogenous   nfac]) ~= 0
                        error('Size mismatch')
                    end
                else
                    if lags > 0
                        warning(['You did not provide a prior mean for the Factor AR coeff. '...
                            'Assume the default values (0.2 on AR1, zeros elsewhere). See the Hitchhiker''s guide.'])
                    end
                end
                % variance
                if isfield(options.priors.F.Phi,'cov') == 1
                    priors.F.Phi.cov   = options.priors.F.Phi.cov;                   
                    if length(priors.F.Phi.cov) ~= (nfac*lags+(1-noconstant) +timetrend+nexogenous ) || size(priors.F.Phi.cov,1 )~=size(priors.F.Phi.cov,2)
                        error('Size mismatch: Covariance Phi should be square, e.g. size(Phi.mean,1) x size(Phi.mean,1)')
                    end
                else
                    if lags > 0
                        warning(['You did not provide a Covariance for the Factor AR coeff. '...
                            'Assume the default values (5 * identity matrix). See the Hitchhiker''s guide.'])
                    end
                end
            else
                if lags > 0                    
                    warning(['You did not provide priors for the AR part of the Factors. '...
                        'Assume the default values. See the Hitchhiker''s guide.'])
                end
            end
        
            % =========================================================== %
            % Priors for the Factor Covariance                        
            % Default values
            %               priors.F.Sigma.scale = eye(nfac); 
            %               priors.F.Sigma.df    = nfac + nexogenous + timetrend + 1;
            if isfield(options.priors.F,'Sigma') == 1
                % scale
                if isfield(options.priors.F.Sigma,'scale') == 1
                    priors.F.Sigma.scale = options.priors.F.Sigma.scale;
                    if size(priors.F.Sigma.scale) ~= [nfac nfac]
                        error('Size mismatch')
                    end
                else
                    if lags > 1 
                        warning(['You did not provide a prior scale for the Factor Covariance. '...
                            'Assume the default values. See the Hitchhiker''s guide.'])
                    end
                end
                % degrees of freedom
                if isfield(options.priors.F.Sigma,'df') == 1
                    priors.F.Sigma.df = options.priors.F.Sigma.df;
                    if length(priors.F.Sigma.df) ~= 1
                        error('Size mismatch')
                    end
                    if priors.F.Sigma.df/2 <= nfac-1
                        error('Too few degrees of freedom - Increase prior df')
                    end
                else
                    if lags > 1 
                        warning(['You did not provide the degrees of freedom for the Factor Covariance. '...
                            'Assume the default values (identity matrix). See the Hitchhiker''s guide.'])
                    end
                end
            else
                if lags > 1                    
                    warning(['You did not provide prior scale and the degrees of freedom for the Factor Covariance. '...
                        'Assume the default values (nfac+1). See the Hitchhiker''s guide.'])
                end
            end
            
            % =========================================================== %
            % Priors for the Loadings
            if isfield(options.priors.F,'Lambda') == 1
                % Default values
                %                 priors.F.Lambda.mean = 0;
                %                 priors.F.Lambda.cov = 10;
                % mean
                if isfield(options.priors.F.Lambda,'mean') == 1
                    priors.F.Lambda.mean  = options.priors.F.Lambda.mean;
                    if length(priors.F.Lambda.mean) ~= 1
                        error('Size mismatch: lambda prior mean should be a scalar; same mean prior for all loadings.')
                    end
                else
                    warning(['You did not provide a prior mean for the factor loadings. '...
                        'Assume the default value (zero). See the Hitchhiker''s guide.'])                    
                end
                % variance
                if isfield(options.priors.F.Lambda,'cov') == 1
                    priors.F.Lambda.cov   = options.priors.F.Lambda.cov;                   
                    if length(priors.F.Lambda.cov) ~= 1
                        error('Size mismatch: covariance Lambda should be a scalar (i.e. the variance)')
                    end
                else
                    warning(['You did not provide prior variance for the factor loadings. '...
                        'Assume the default values (10). See the Hitchhiker''s guide.'])   
                end
            else
                warning(['You did not provide priors for the factor loadings. '...
                    'Assume the default values (zero mean and 10 variance). See the Hitchhiker''s guide.'])               
            end             
        else
            warning(['You did not provide any prior for the factors (loadings, AR and Sigma). '...
                'Assume the default values. See the Hitchhiker''s guide.'])                                    
        end    
        % =========================================================== %
        % Priors for the idiosyncratic term
        % =========================================================== %
        if isfield(options.priors,'G') == 1
            % Default Values
            %             priors.G.Sigma.scale = 1;
            %             priors.G.Sigma.df    = 4;
            % =========================================================== %
            % priors for the standard devitiations
            if isfield(options.priors.G,'Sigma') == 1
                % scale
                if isfield(options.priors.G.Sigma,'scale') == 1
                    priors.G.Sigma.scale = options.priors.G.Sigma.scale;
                    if length(priors.G.Sigma.scale) ~= 1
                        error('Size mismatch')
                    end
                else
                    warning(['You did not provide a prior scale for the idiosyncratic error standard deviations.'...
                        'Assume the default value (1). See the Hitchhiker''s guide.'])
                end
                % degrees of freedom
                if isfield(options.priors.G.Sigma,'df') == 1
                    priors.G.Sigma.df = options.priors.G.Sigma.df;
                    if length(priors.G.Sigma.df) ~= 1
                        error('Size mismatch')
                    end
                    if priors.G.Sigma.df/2 <= 1
                        error('Too few degrees of freedom - Increase prior df')
                    end
                else
                    warning(['You did not provide the degrees of freedom for the idiosyncratic error standard deviations.'...
                        'Assume the default values (4). See the Hitchhiker''s guide.'])
                end
                
            else
                warning(['You did not provide prior scale and the degrees of freedom for the idiosyncratic error standard deviations.'...
                    'Assume the default values (1 and 4). See the Hitchhiker''s guide.'])
            end
            if isfield(options.priors.G,'Phi') == 1
                warning(['AR1 for the idiosyncratic errors is currently not supported by the toolbox.'...
                    'Assume that all AR1 coefficients of the idiosyncratic errors are zero.'])
            end
%             if ilags == 0
%                 error('If you specify the priors for the idiosyncratic error (options.priors.G), you need also to specify the number of lags (options.ilags)')
%             end
        end
        %======================================================================
        % Jeffrey priors on Factor dynamics
        %======================================================================
        if isfield(options.priors,'name')==1 && strmatch(options.priors.name,'Jeffrey') == 1
            dummy = 0;
            priors.name= 'Conjugate-and-Jeffry-on-Factors-Dynamics';
        end

    end


    %======================================================================
    % IRF options
    %======================================================================
    if isfield(options,'hor') ==1
        hor = options.hor;
    end
    if isfield(options,'long_run_irf')==1
        % Activating Long run IRF
        long_run_irf = options.long_run_irf;
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
        narrative_signs_irf = 1; 
        narrative           = options.narrative ;
    end
    if isfield(options,'proxy')==1
        % Activating IRF with provy
        proxy_irf   = 1;
        in.proxies  = options.proxy;
        %in.vars     = y;
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
    end

end

%********************************************************
%* Prepare the data
%********************************************************
if do_the_demean == 1
    y = demean(y);
end

%********************************************************
%* Consistency Checks
%********************************************************
if lags < l_F
    error('The number of lags in the factor transition equation must be as large as the number of lagged factor entering the measurement equation');
end

% constant in the factor equations
nx    = 1-noconstant; 

% number of regressors in the factor equation
nk    = nfac*lags + nx + timetrend + nexogenous; 

% number of factors
params.K_F = nfac;

% number of variables/idiosyncratic errors
params.K_G = ny;

% lag order for factors and idiosyncratic errors
params.q_F = lags;
params.q_G = ilags;

% factors/idio entering with lags in the measurement equation (default 0)
params.l_F = l_F;
params.l_G = l_G;


%********************************************************
%* Preliminary PC analysis to initialize the Gibbs sampler
%********************************************************

% paramsF.B        = 1; 
paramsF.q_G      = ilags;   % lag structure of idiosyncratic errors
paramsF.K_blocks = ny;      % number of idiosyncratic errors
paramsF.q_F      = lags*ones(nfac,1); 
paramsF.K_facs   = nfac; 
paramsF.l_F      = l_F;

% Initialize factors as PC estimates plus random noise
[~,f,~,~] = pc_T(y,nfac,0);
f         = f + noise_var*randn(T,nfac);
for k = 1 : nfac
    f(:,k) = sign(corr(f(:,k),y(:,k)))*f(:,k);
end 

% Now initialize parameters based on starting values for the factors
% First draw starting values of aggregate factors on their own lags
% to obtain starting values for the psi_F's (autoregressive factor) and
% sig2_F's (variance of factor) 
% Draw VAR-COV parameters of the factors (if dynamic)
if lags > 0 
    [psi_F, sig2_F] = draw_Psi_Sigma(f, priors.F, lags);
end

% Now regress starting values of block factors on starting values of
% aggregate factors to obtain starting value for Lambda_F
F_lags = f;
for ell = 1 : max(l_F)
    %F_lags = [F_lags lagn(f,ell)];
    F_lags = [F_lags lagX(f,ell)];
end
lambda_F = ((F_lags'*F_lags)\(F_lags'*y))';
for ell = 1 : max(l_F)+1
    Lambda_F{ell} = lambda_F(:, (ell-1)*nfac+1 : ell*nfac);        
end
for i1 = 1 : nfac                
    signlam(i1) = sign(Lambda_F{1}(i1,i1));
    f(:,i1) = f(:,i1)*signlam(i1);
    for ell = 1 : l_F+1
        Lambda_F{ell}(:,i1) = Lambda_F{ell}(:,i1)*signlam(i1);
    end
end     
% Fix upper left block of Lambda_F{1} to be lower triangular with ones on diagonal
for j1 = 1 : nfac
    Lambda_F{1}(j1,j1) = 1;
    Lambda_F{1}(j1,j1+1:nfac) = zeros(1,nfac - j1);
end

% lambda_Fb       = cell2mat(Lambda_F);
Lambda_F_mat    = cell2mat(Lambda_F);
eG              = compute_resids(y, f, Lambda_F_mat, paramsF); % oo = y - f* Lambda_F_mat';

% paramsGb.q_F    = ilags*ones(ny,1); 
% paramsGb.K_facs = ny; 
% paramsGb.l_F    = l_G;
% draw starting values of idiosyncratic part on their own lags
% to obtain starting values for the psi_G's (autoregressive factor) and
% sig2_G's (variance of factor) 
% start from diffuse priors
if ilags == 0
    psi_G  = zeros(ny,1);
    sig2_G = diag(eG'*eG./T);
else
    psi_G  = zeros(ny,1);
    sig2_G = ones(ny,1);
    for ww = 1 : ny
        [psi_G(ww,1), sig2_G(ww,1)] = draw_Psi_Sigma(eG(:,ww), priors.G, ilags); 
    end
    %[psi_G,sig2_G] = draw_rho_sigma(eG,zeros(ny,1),0.5*ones(ny,1),startprior,paramsGb,sig_fix);
end

% Settings for the forecasts
forecast_data.xdata       = ones(fhor, nx);
if timetrend
    forecast_data.xdata = [forecast_data.xdata (T-lags+1 : T-lags+fhor)'];
end

%**************************************************
%* Generating draws form the Posterior Distribution
%**************************************************
% Preallocation of memory
% Matrices for collecting draws from Posterior Density
% Last dimension corresponds to a specific draw

f_draws      = zeros(T,nfac,K);         % factors (smoothed)
f_filt_draws = zeros(T,nfac,K);         % factors (filtered)
ef_draws     = zeros(T-lags,nfac,K);    % factor innovations
eg_draws     = zeros(T,ny,K);           % idiosyncratic errors
lambda_draws = zeros(ny,nfac,K);        % loadings
Phi_draws    = zeros(nk,nfac,K);        % factor VAR Autoregressive coefficient
Sigma_draws  = zeros(nfac,nfac,K);      % factor VAR Covariance 
phi_draws    = zeros(ny,K);             % idiosyncratic errors AR
sigma_draws  = zeros(ny,K);             % idiosyncratic errors Variance
ir_draws     = zeros(ny,hor,nfac,K);    % variable, horizon, shock and draws - Cholesky IRF  
irlr_draws   = zeros(ny,hor,nfac,K);    % variable, horizon, shock and draws - Long Run IRF
Qlr_draws    = zeros(nfac,nfac,K);      % long run impact matrix
if signs_irf == 1
    irsign_draws = ir_draws;
    Omega_draws  = Sigma_draws;
end
if narrative_signs_irf == 1
    irnarrsign_draws = ir_draws;
    Omegan_draws     = Sigma_draws;
end
if proxy_irf == 1
    irproxy_draws = ir_draws;
end
% if zeros_signs_irf == 1
%     irzerosign_draws   = ir_draws;
%     Omegaz_draws       = Sigma_draws;
% end
yhatfut_no_shocks         = NaN(fhor, ny, K);   % forecasts with factor-shocks (idiosyncratic shocks are zero)
yhatfut_with_shocks       = NaN(fhor, ny, K);   % forecast without the factor-shocks (idiosyncratic shocks are zero)

%**************************************************
%* Gibbs Sampler 
%**************************************************
waitbar_yes = 0;
if K > 99
    waitbar_yes = 1;
    wb = waitbar(0, 'Generating draws from the Posterior Distribution');
end
do_pagemtimes = 0;
if exist('pagemtimes','builtin') == 5
    do_pagemtimes = 1;
end
d = 0; 
for d1 = 1 : skip_*K + burnin_
    %======================================================================
    % Inferece: Drawing from the posterior distribution
    % 1: draw factors 
    if lags == 0 
        % static factor model
        f_filt      = zeros(T,nfac);
        iSigmay     = (Lambda_F_mat * Lambda_F_mat' + diag(sig2_G))\eye(ny);        %iSigmay = inv(Lambda_F_mat * Lambda_F_mat' + diag(sig2_G));
        f_var       = eye(nfac) -  Lambda_F_mat' * iSigmay * Lambda_F_mat;
        f_var_chol  = chol(f_var,'lower');
        i_f         = randn(T,nfac);
        f_mean      = y * iSigmay' * Lambda_F_mat;
        f           = f_mean + i_f * f_var_chol';                            
    else
        % dynamic factors (Kalman Smoother)
        %f              = sample_facs(y, zeros(size(f)), Lambda_F_mat, psi_F, sig2_F, psi_G, sig2_G, paramsF);
        %[f, f_filt]    = sample_facs(y, zeros(size(f)), Lambda_F_mat, psi_F, sig2_F, psi_G, sig2_G);
        [f , f_filt]   = sample_facs(y, zeros(size(f)), Lambda_F_mat, psi_F, sig2_F, psi_G, sig2_G);
    end
    
    % 2: draw loadings 
    Lambda_F_mat   = draw_lambda(y, f, psi_G, sig2_G, priors.F.Lambda, paramsF);    
    
    % 3: draw VAR-COV parameters of the factors
    if lags == 0 % static factor model
        psi_F  = zeros(nk,nfac);
        sig2_F = eye(nfac);
    else % dynamic factors (MN-IW)
        [psi_F, sig2_F] = draw_Psi_Sigma(f, priors.F, lags);
    end
    
    % 4: draw idiosyncratic errors
    eG              = compute_resids(y, f, Lambda_F_mat, paramsF); 
    
    % 5: draw ar/var parameters of the idiosyncratic errors
    if ilags ~= 0        
        for ww = 1 : ny
            [psi_G(ww,1), sig2_G(ww,1)] = draw_Psi_Sigma(eG(:,ww), priors.G, ilags);
        end
%         %[psi_G,sig2_G]  = draw_rho_sigma(eG, psi_G, sig2_G, priors.G, paramsGb, sig_fix);
     else
%         sig2_G = diag(eG'*eG./T);
          sig2_G = draw_sigma(eG, f, Lambda_F_mat, priors.G);
    end
    % store the draws
    if d1 > burnin_ && mod(d1,skip_) == 0 
        d = d+1;
        f_draws(:,:,d)      = f;
        f_filt_draws(:,:,d) = f_filt;
        eg_draws(:,:,d)     = eG;
        lambda_draws(:,:,d) = Lambda_F_mat;
        Phi_draws(:,:,d)    = psi_F;
        Sigma_draws(:,:,d)  = sig2_F;
        phi_draws(:,d)      = psi_G;
        sigma_draws(:,d)    = sig2_G;
        [ff,FF]             = YXB_(f,lags,[nx timetrend]);
        errors              = ff - FF * psi_F;
        ef_draws(:,:,d)     = errors;

        if lags > 0
            % IRF and Forecasts are computed only with dynamic factors
            %======================================================================
            % IRF
            % Compute the impulse response functions
            %C_       = rescaleFAVAR(STD,Lambda,size(y1,2),order_pc);
            % with cholesky
            tmp                    = iresponse(psi_F,sig2_F,hor,eye(nfac));
            if do_pagemtimes == 1
                Lam_                   = repmat(Lambda_F_mat,1,1,nfac);
                ir_draws(:,:,:,d)      = pagemtimes(Lam_,tmp);
            else
                for ff = 1 : nfac
                    ir_draws(:,:,ff,d) = Lambda_F_mat * tmp(:,:,ff);
                end
            end
            % with long run restrictions
            if long_run_irf == 1
                [irlr,Omega]           = iresponse_longrun(psi_F,sig2_F,hor,lags);
                if do_pagemtimes == 1
                    irlr_draws(:,:,:,d)    = pagemtimes(Lam_,irlr);
                else
                    for ff = 1 : nfac
                        irlr_draws(:,:,ff,d) = Lambda_F_mat * tmp(:,:,ff);
                    end
                end
                Qlr_draws(:,:,d)       = Omega;
            end
            % with sign restrictions
            if signs_irf == 1
                [irsign,Omega]         = iresponse_sign(psi_F,sig2_F,hor,signs,Lambda_F_mat);
                irsign_draws(:,:,:,d)  = irsign;
                Omega_draws(:,:,d)     = Omega;
            end
            % with narrative and sign restrictions
            if narrative_signs_irf == 1
                [irnarrsign,Omega]         = iresponse_sign_narrative(errors,psi_F,sig2_F,hor,signs,narrative,Lambda_F_mat);
                irnarrsign_draws(:,:,:,d)  = irnarrsign;
                Omegan_draws(:,:,d)        = Omega;
            end
            % with proxy
            if proxy_irf == 1
                in.Phi                  = psi_F;
                in.Sigma                = sig2_F;
                for nf = 1 : nfac
                    nfnot                   = setdiff(1 : nfac, nf);
                    in.res                  = ef_draws(:, [nf, nfnot], d);
                    in.vars                 = f;
                    tmp_                    = iresponse_proxy(in);
                    irproxy_draws(:,:,nf,d) = Lambda_F_mat * tmp_.irs';
                    clear tmp_
                end
            end
            
            %======================================================================
            % Forecasts
            % compute the out of sample forecast (unconditional)
            forecast_data.initval     = f(end-lags+1:end, :);
            [frcst_no_shock,frcsts_with_shocks] = forecasts(forecast_data,psi_F,sig2_F,fhor,lags);
            yhatfut_no_shocks(:,:,d)            = frcst_no_shock * Lambda_F_mat';
            yhatfut_with_shocks(:,:,d)          = frcsts_with_shocks * Lambda_F_mat';
        end
    end
    if waitbar_yes, waitbar(d1/(skip_*K  + burnin_), wb); end
end
if waitbar_yes, close(wb); end

%==========================================================================
% Collecting the output

BDFM.lags       = lags;               % lags factors
BDFM.ilags      = ilags;              % lags idio
BDFM.N          = ny;                 % number of variables
BDFM.nfac       = nfac;               % number of factors
BDFM.prior      = priors;             % priors used
BDFM.params     = params;             %
BDFM.varnames   = varnames;           % variables names
BDFM.ndraws     = K;

% Inference: Reduced Form Parameters
% the last dimension of these objects corresponds to a draw from the posterior
BDFM.Phi_draws    = Phi_draws;          % draws from the autoregressive part (Factor)
BDFM.Sigma_draws  = Sigma_draws;        % draws from the covarance matrix (Factor)
BDFM.lambda_draws = lambda_draws;       % draws from the factor loading matrix
BDFM.phi_draws    = phi_draws;          % draws from the AR coeff of the idiosyncratic shocks
BDFM.sigma_draws  = sigma_draws;        % draws from the variance of the idiosyncratic shocks 
BDFM.f_draws      = f_draws;            % draws from the factors (smoothed)
BDFM.f_filt_draws = f_filt_draws;       % draws from the factors (filtered)

BDFM.eg_draws     = eg_draws;           % draws from the idiosyncratic errors
BDFM.ef_draws     = ef_draws;           % draws from the innovation in the factors

% Inference: prediction (forecast horizon, variable)
BDFM.fhor                     = fhor;                      % forecast horizon
BDFM.forecasts.no_shocks      = yhatfut_no_shocks;         % trajectories of forecasts without shocks
BDFM.forecasts.with_shocks    = yhatfut_with_shocks;       % trajectories of forecasts with shocks

% Inference: IRFs (variable, horizon, shock and draws)
BDFM.hor          = hor;                % IRF horizon
BDFM.ir_draws     = ir_draws;           % draws from the IRF with cholesky
BDFM.irlr_draws   = irlr_draws;         % draws from the IRF with Long Run
BDFM.Qlr_draws    = Qlr_draws;          % Long Run Rotation matrix
if signs_irf == 1 && narrative_signs_irf == 0
    BDFM.irsign_draws = irsign_draws;
    BDFM.Omegas       = Omega_draws;
else
    BDFM.irsign_draws = [];
    BDFM.Omegas       = [];
end
if narrative_signs_irf == 1
    BDFM.irnarrsign_draws = irnarrsign_draws;
    BDFM.Omegan           = Omegan_draws;
else
    BDFM.irnarrsign_draws = [];
    BDFM.Omegan           = [];
end
if proxy_irf == 1
    BDFM.irproxy_draws = irproxy_draws;
else
    BDFM.irproxy_draws= [];
end
% if zeros_signs_irf == 1
%     BDFM.irzerosign_draws   = irzerosign_draws;
%     BDFM.Omegaz             = Omegaz_draws;
% else
%     BDFM.irzerosign_draws = [];
%     BDFM.Omegaz           = [];
% end


%% Child Functions

% child function 1: Draw AR and STD from normal-inverse gamma (independent processes)
function [PSI,SIG] = draw_rho_sigma(Ft,psi_F,sig2_F,prior,params,var_fix)

    q_F     = params.q_F; %
    max_q_F = max(q_F);
    max_l_F = max(l_F);
    PSI     = psi_F;
    SIG     = sig2_F;
    K_F     = size(Ft,2);
    
    for k1 = 1 : K_F
        q_Fk    = q_F(k1);
        xpsi_F1 = zeros(max_q_F,1);
        if q_Fk > 0
            ystar = Ft(:,k1);
            xstar = []; xstar1 = [];
            for i= 1 : q_Fk
                xstar = [xstar lagX(Ft(:,k1),i)];             
                %xstar  = [xstar lagn(Ft(:,k),i)];
                %xstar1 = [xstar1 lagmatrix(Ft(:,k),i)];
                %test   = max(max(abs(xstar-xstar1)));
            end
            ystar = trimr(ystar,max(max_q_F,max_l_F),0);
            xstar = trimr(xstar,max(max_q_F,max_l_F),0);
            TT    = size(ystar,1);
            
            % bring to right dimension
            prior_mean    = prior.Phi.mean*ones(q_Fk,1);
            inv_prior_var = inv(prior.Phi.cov)*eye(q_Fk);
            sig_inv = 1/sig2_F(k1);
            
            % posterior mean and variance
            post_var = inv(inv_prior_var + sig_inv*xstar'*xstar);
            post_mean = post_var*(inv_prior_var*prior_mean + sig_inv*xstar'*ystar);
            
            % draw from multivariate normal with posterior mean and variance
            C = chol(sig2_F(k1)*post_var);
            accept = 0; counter = 0;
            while accept ==0 && counter < 1000
                psi_F1  = post_mean+C'*randn(q_Fk,1);
                ceof    = flipud([1; -psi_F1]);
                root    = roots(ceof);
                rootmod = abs(root);
                if min(rootmod) > 1.0001
                    accept=1;
                else
                    counter=counter+1;
                    accept=0;
                end
            end % while
            if ~accept && counter >= 1000
                disp(counter);
                return;
            end
            SSE = (ystar-xstar*psi_F1)'*(ystar-xstar*psi_F1);
            xpsi_F1(1:q_Fk) = psi_F1(1:q_Fk);
        else % if q_Fk >0
            
            TT    = size(Ft,1); %TT    = rows(Ft);
            ystar = Ft(:,k1);
            SSE   = ystar'*ystar;
            
        end % q_Fk >00
        
        if ~isempty(var_fix)
            sig2_F1 = var_fix;
        else
            dd       = prior.Sigma.scale + SSE;
            c       = chi2rnd(TT + prior.Sigma.df,1);
            sig2_F1 = dd/c;
        end % if DO_FIX_VARIANCE
        
        %xpsi_F1(1:q_Fk) = psi_F1(1:q_Fk);
        if k1 == 1
            SIG = sig2_F1;
            PSI = xpsi_F1';
        else
            SIG = [SIG; sig2_F1];
            PSI = [PSI; xpsi_F1'];
        end        
    end % for k = 1:K_F
end % end child function 1

% child function 2: Draw Phi and Sigma from MultivariateNormal-InverseWishart
function [Phi,Sigma] = draw_Psi_Sigma(yfac,prior,lags)
    
        ydata = yfac;
        %xdata = ones(T,nfac);
        xdata = [];
        %********************************************************
        % possible future options to be added
        %********************************************************
        %if timetrend ==1
        %    xdata = [xdata [1-lags : T-lags]'];
        %    % xdata = [xdata [1:T]'];
        %end
        %if nexogenous >0
        %    xdata = [xdata exogenous(idx,:)];
        %end        
        %if heterosked == 0
        %    var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu);
        %else
        %    var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu, ww);
        %end
        %********************************************************
        % OLS 
        %********************************************************
        var = rfvar3(ydata, lags, xdata, [T; T], 0, 0);
        Tu  = size(var.u, 1);
        Ny  = size(var.y, 2);
        Nx  = size(var.X, 2);
        
        if dummy == 0
            %********************************************************
            % Flat Jeffrey Prior
            %********************************************************
            flat            = 1;
            posterior.df    = Tu - ny*lags - nx + flat*(ny+1) - nexogenous;
            posterior.S     = var.u' * var.u;
            posterior.XXi   = var.xxi;
            posterior.PhiHat = var.B;
        
        else
            %********************************************************
            % Step 0: Compute the Posterior Distribution Moments (MN-IW) 
            %********************************************************
            Ai              = inv(prior.Phi.cov);
            %posterior.df    = Tu - Ny*lags - nx - nexogenous - prior.Sigma.df;
            posterior.df    = Tu - Nx - prior.Sigma.df;
            posterior.XXi   = inv(var.X'*var.X + Ai);
            posterior.PhiHat = posterior.XXi * (var.X' * var.y + Ai * prior.Phi.mean);
            posterior.S = var.u' * var.u + prior.Sigma.scale + ...
                prior.Phi.mean' * Ai * prior.Phi.mean + ...
                var.B' * (var.X'*var.X) * var.B ...
                - posterior.PhiHat' * (var.X'*var.X + Ai) * posterior.PhiHat;
        end
        %********************************************************
        % Inferece: Draw Phi-Sigma from the posterior distribution
        %********************************************************
        % Step 1: Draw from the Covariance
        S_inv_upper_chol = chol(inv(posterior.S));
        Sigma            = rand_inverse_wishart(Ny, posterior.df, S_inv_upper_chol);
        %********************************************************
        % Step 2: given the Covariance Matrix, draw from the AR parameters
        Sigma_lower_chol = chol(Sigma)';
        XXi_lower_chol   = chol(posterior.XXi)';
        Phi1 = randn(Nx * Ny, 1);
        Phi2 = kron(Sigma_lower_chol , XXi_lower_chol) * Phi1;
        Phi3 = reshape(Phi2, Nx, Ny);
        Phi  = Phi3 + posterior.PhiHat;        
        
end

% child function 3: compute idiosyncratic residuals
function e = compute_resids(obs,fac,Lambda,params)
    
    K_F     = params.K_facs;
    % l_F     = params.l_F;
    max_l_F = max(l_F);
    Fmat    = fac;
    for j0 =1 : max_l_F
        %Fmat = [Fmat lagn(fac,j0)];
        Fmat = [Fmat lagX(fac,j0)];
    end
    fit = zeros(size(obs));
    for j0 = 1 : K_F
        indx  = j0 : K_F : K_F*(max_l_F+1);
        fit   = fit + Fmat(:,indx)*Lambda(:,indx)';
    end
    e = obs-fit;
end % end child funtion 3

% child function 4: draw the unobserved factors 
function [FT_mat, f_filt] = sample_facs(g,alpha,lambda_tilde,psi_F,sig2_F,psi_G,sig2_G)
  
    Gmat       = g;
    sig2_G_vec = sig2_G;
    psi_G_mat  = psi_G;
    
    if any(psi_G_mat ~= 0) 
        ystar = zeros(T,ny);
        for i = 1 : ny
            ystar(:,i) = filter([1 -psi_G_mat(i,:)],1,Gmat(:,i));
        end
        ystar = ystar';
    else    
        ystar = Gmat';
    end
    
    % transition matrix of the filtered data is r= p+1+q+1    
    % companion form 
    A       = [psi_F(1 : nfac * lags, :)'; eye(nfac*(lags-1), nfac*lags)];
    H       = [lambda_tilde, zeros(size(lambda_tilde,1) , nfac*(lags-1-l_F))];
    ndim    = size(H,2);
    GG      = eye(nfac * lags, nfac);    
    QQ      = GG * sig2_F * GG';
    R       = diag(sig2_G_vec);
    Alpha   = [alpha zeros(T , nfac*(lags-1))];
    
    % initialize Kalman filter
    Gtt     = zeros(nfac * lags, 1);
    Ptt     = lyapunov_symm(A, QQ);
%     Ptt     = P00;%(1:nfac , 1:nfac);
%     P00     = inv(eye(ndim^2)-kron(A,A))*vec(sig2_F);
%     Ptt     = reshape(P00,ndim,ndim);
    Fmat  = zeros(nfac * lags,T);
    Pmat  = zeros(nfac * lags,nfac * lags,T);
    
    if use_the_ss_kalman_gain == 1
        % compute the SS kalman gain
        [SigSS,KapSS] = GetSigKap(A,H,QQ,R);
    end
    
    % Kalman filter recursion
    for t = 1 : T
        Gtt1   = Alpha(t,:)' + A*Gtt;         % x(t|t-1)
        Ptt1   = A*Ptt*A' + QQ;               % P(t|t-1)
        ett1   = ystar(:,t) - H*Gtt1;         % eta(t|t-1) = y(t)- y(t|t-1)= y(t)-H*x(t|t-1)
        if use_the_ss_kalman_gain
            k_gn   = KapSS;
            Ptt    = SigSS;
        else
            v_ett1 = H*Ptt1*H' + R;           % var(eta(t|t-1))
            k_gn   = Ptt1*H'/v_ett1;          % K      = P(t|t-1)H'/ f(t|t-1)
            Ptt    = (eye(ndim)-k_gn*H)*Ptt1;     % P(t|t) = P(t|t-1)-K*H*P(t|t-1)
        end
        Gtt         = Gtt1 + k_gn*ett1;            % x(t|t) = x(t|t-1)+ K eta(t|t-1)        
        Fmat(:,t)   = Gtt;       
        Pmat(:,:,t) = Ptt;                    % Pmat{t}     = Ptt;
    end
    f_filt = Fmat(1:nfac,:)';

    % Carter-Kohn backward sampling algorithm
    FT_mat          = zeros(T, nfac);         
    %G               = chol(Pmat{T})';
    [G,flag]        = chol(Pmat(:,:,T),'lower');
    if flag ~= 0, keyboard; end        
    FT              = Fmat(:,T)+G*randn(ndim,1);
    FT_mat(T, 1:nfac) = FT(1:nfac)';           
    for t= T-1 : -1 : 1
        jj              = 1 : nfac;          
        etT             = FT(jj) - Alpha(t+1,jj)' - A(jj,:)*Fmat(:,t);
        v_etT           = A(jj,:)* Pmat(:,:,t)*A(jj,:)' + QQ(jj,jj);
        k_gn0           = Pmat(:,:,t) * A(jj,:)' /v_etT;
        Fmat(:,t)       = Fmat(:,t) + k_gn0*etT;
        Pmat(:,:,t)     = (eye(ndim)-k_gn0*A(jj,:))*Pmat(:,:,t);
        G               = chol(Pmat(:,:,t))';
        FT              = Fmat(:,t) + G*randn(ndim,1);
        FT_mat(t, jj)   = FT(jj, 1)';
    end % end t   
    
end % end child function 4

% child function 5: draw the factor loadings
function LAMBDA_mat = draw_lambda(G,F,psi_G,SIG2_G,prior,params)

    q_G     = params.q_G;
    max_q_G = max(q_G);
    % K_G     = params.K_blocks;
    K_F     = params.K_facs;
    l_F     = params.l_F;
    max_l_F = max(l_F);
    Lamda_F = [];
    count   = 0;
    Gb      = G;
    [~,k_G] = size(Gb);
    l_Fb    = l_F;
    
    for i= 1 : k_G
        
        count = count + 1;
        F_filt = filter([1 -psi_G(i,:)], 1, F);
        ystar  = filter([1 -psi_G(i,:)], 1, Gb(:,i));
        Fstar  = F_filt; %Fstar1  = F_filt;
        for j = 1 : l_Fb
            %Fstar = [Fstar lagn(F_filt,j)];
            Fstar = [Fstar lagX(F_filt,j)];
            %tolll = max(max(abs(Fstar1-Fstar)));
        end
        Fstar = trimr(Fstar, max(max_q_G,l_Fb), 0);
        ystar = trimr(ystar, max(max_q_G,l_Fb), 0);
        
        % bring to right dimension
        prior_mean = prior.mean*ones(K_F*(l_Fb+1),1);
        prior_var  = prior.cov*eye(K_F*(l_Fb+1));
        sig_inv    = 1/SIG2_G(i);
        
        % posterior mean and variance
        post_var  = inv(inv(prior_var) + sig_inv*Fstar'*Fstar);
        post_mean = post_var*(inv(prior_var)*prior_mean + sig_inv*Fstar'*ystar);
        
        % draw from multivariate normal with posterior mean and variance
        C        = chol(post_var);
        lamda_Fi = post_mean + C'*randn(K_F*(l_Fb+1),1);
        Lamda_F  =[Lamda_F; lamda_Fi' zeros(1,max_l_F-l_Fb)];
        
    end % end i
    
    % Impose lower triangularity with ones on the diagonal for identification
    for l = 1 : max_l_F+1
        LAMBDA{l} = Lamda_F(:,(l-1)*K_F+1:l*K_F);
        for j = 1 : K_F
            if l == 1
                LAMBDA{l}(j,j) = 1;
            end
            LAMBDA{l}(j,j+1:K_F) = zeros(1,K_F-j); % make loadings lower triangular at all lags
        end
    end
    
     LAMBDA_mat   = cell2mat(LAMBDA);    
    
end % end child function 5

% child function 6: matrix manipulation
function xx=trimr(x,a,b)
[nt,~] = size(x);
if a >= 0 
    xx=x(a+1:nt,:);
end
if b >=0 
    if a > 0 
        x=xx; 
    end
    [nt,nc]=size(x);
    xx=x(1:nt-b,1:nc);
end
end% end child function 6

% child function 7: Draw STD from inverse gamma (independent processes)
function sig2_G = draw_sigma(eG, f, LAMBDA, priors)
    
    sig2_G      = nan(ny,1);
    L_var_prior = eye(nfac);
%     L_OLS       = (f'*f)\(f'*y(:,nfac+1 : ny));
%     %R_OLS = (y - f*LAMBDA')'*(X - f*LAMBDA')./(T-ny);        
%     L=[eye(nfac) L_OLS]';
    L = LAMBDA;
    alpha=0.001;
    %s0=3;
    s0 = priors.Sigma.scale;
    for n=1 : ny        
        % draw sigma_G
        CC    = L_var_prior+(f'*f)\eye(size(f'*f));
        R_bar = s0 + eG(:,n)'*eG(:,n) + L(n,:)*(CC\L(n,:)');       % R_bar=s0+ed'*ed+L(n,:)*inv(L_var_prior+inv(FY'*FY))*L(n,:)';
        Rd    = chi2rnd(T+alpha);
        Rd=R_bar/Rd;          %R(n,n)=Rd;
        sig2_G(n) = Rd;
    end


end

end % end bdfm_ 
