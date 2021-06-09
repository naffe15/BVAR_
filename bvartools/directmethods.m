function [dm] = directmethods(y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'directmethods' generates local projeciton IRF and direct forecasts

% Core Inputs:
% - y, data columns variables
% - lags, lag order of the VAR

% Additonal Inputs collected options: see below

% Filippo Ferroni, 27/02/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('the BVAR toolbox needs at least two inputs: data and number of lags');
end
if lags < 1
    error('lags cannot be zero or negative');
end
% number of observable variables
[T,ny]       = size(y);

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
     pkg load statistics
    randn('state',999);
    rand('state',999);
end

hor          = 24;
conf_sig     = 0.9;
controls_    = 0;
robust_se_   = 1; % by default robust SE
proxy_       = 0;
noconstant   = 0;
ns           = 0;
Q            = eye(ny);
dummy        = 0;
K            = 5000;
max_prior_tau_ = 0;
lb             = -1e10;
ub             = 1e10;


%********************************************************
%* CUSTOMIZED SETTINGS
%********************************************************
if nargin>2
    %======================================================================
    % Various options
    %======================================================================
    if isfield(options,'hor') ==1 % horizon of IRF and Forecasts
        hor = options.hor;
    end
    if isfield(options,'conf_sig') ==1 % CI of OLS estimation of LP and DF
        conf_sig = options.conf_sig;
    end
    if isfield(options,'Q') == 1 % orthonormal matrix
        Q = options.Q;
    end
    if isfield(options,'K') == 1 % # of draws for the BLP
        K = options.K;
    end
    %======================================================================
    % options: adding controls
    %======================================================================
    if isfield(options,'controls') ==1
        controls_  = 1;
        controls   = options.controls;
        if T~=size(controls,1)
            error('Control variables and observables must have the same time length')
        end
    end
    %======================================================================
    % options: adding proxy variable for identification
    %======================================================================
    if isfield(options,'proxy') ==1
        proxy_  = 1;
        proxy   = options.proxy;
        ns      = size(proxy,2);
        if T~=size(proxy,1)
            error('Shocks proxies and observables must have the same time length')
        end
    end
    %======================================================================
    % options: Compuing Robust Standard Errors
    %======================================================================
    if isfield(options,'robust_se_') == 1
        robust_se_ = options.robust_se_;
        % robust_se_ = 0    unadjusted SE
        % robust_se_ = 1    NW Robust SE:  Hamilton (1994), Ch 10 pag 282, eq (10.5.20)
        % robust_se_ = 5    Matlab HAC function: need Matlab Econ Toolbox
    end
    %======================================================================
    % Conjugate/Hierachical prior options
    %======================================================================
    if (isfield(options,'priors')==1 && strcmp(options.priors.name,'Conjugate')==1) || (isfield(options,'priors')==1 && strcmp(options.priors.name,'conjugate')==1) || ...
            (isfield(options,'prior')==1 && strcmp(options.prior.name,'Conjugate')==1) || (isfield(options,'prior')==1 && strcmp(options.prior.name,'conjugate')==1)
        
        dummy = 2;
        prior.name= 'Conjugate';
        % Priors for the AR parameters
        if isfield(options.priors,'Phi') == 1
            % mean
            if isfield(options.priors.Phi,'mean') == 1
                prior.Phi.mean  = options.priors.Phi.mean;
                if size(prior.Phi.mean) ~= [ny*lags+(1-noconstant)    ny]
                    error('Size mismatch')
                end
            else
                warning(['You did not provide a prior mean for the AR coeff. ' ...
                    'Assume zeros everywhere.'])
                prior.Phi.mean  = zeros(ny*lags+(1-noconstant) , ny);
            end
            % variance
            if isfield(options.priors.Phi,'cov') == 1
                prior.Phi.cov   = options.priors.Phi.cov;
                if length(prior.Phi.cov) ~= (ny*lags+(1-noconstant) ) || size(prior.Phi.cov,1 )~=size(prior.Phi.cov,2)
                    error('Size mismatch: Covariance Phi should be square, e.g. size(Phi.mean,1)x size(Phi.mean,1)x')
                end
            else
                warning(['You did not provide a Covariance for the AR coeff. ' ...
                    'Assume 10 times Identity Matrix.'])
                prior.Phi.cov  = 10 * eye((ny*lags+(1-noconstant) ));
            end
        else
            warning(['You did not provide prior mean and covariance for the AR coeff ' ...
                'Assume zeros everywhere with covariance 10 times Identity Matrix.'])
            %             prior.Phi.cov   = 10 * eye((ny*lags+(1-noconstant) +timetrend ) * ny);
            prior.Phi.cov   = 10 * eye((ny*lags+(1-noconstant)  ));
            prior.Phi.mean  = zeros(ny*lags+(1-noconstant) ,  ny);
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
                warning(['You did not provide a prior mean for the Residual Covariance. ' ...
                    'Assume identity matrix.'])
                prior.Sigma.scale = eye(ny);
            end
            % degrees of freedom
            if isfield(options.priors.Sigma,'df') == 1
                prior.Sigma.df = options.priors.Sigma.df;
                if length(prior.Sigma.df) ~= 1
                    error('Size mismatch')
                end
                if prior.Sigma.df < ny
                    error('Too few degrees of freedom - prior on variance residuals')
                end
            else
                warning(['You did not provide the degrees of freedom for the Residual Covariance. ' ...
                    'Assume ny+1.'])
                prior.Sigma.df = ny + 1;
            end
        else
            warning(['You did not provide prior mean and variance for the Residual Covariance. ' ...
                'Assume an identity matrix matrix with N+1 degrees of freedom.'])
            prior.Sigma.scale = eye(ny);
            prior.Sigma.df    = ny + 1;
        end
        if isfield(options.priors,'tau') == 1 && isnumeric(options.priors.tau) == 1
            prior.tau = options.priors.tau;
            if length(options.priors.tau) ~= hor
                error('tau must be a vector of size hor')
            end
        elseif isfield(options.priors,'max_tau') == 1 && options.priors.max_tau ==1
            max_prior_tau_ = 1;
            max_compute    = 3;
            prior.tau = ones(hor,1);
        else
            prior.tau = ones(hor,1);
        end
    end
    % options for the maximization
    if isfield(options,'max_compute') == 1
        max_compute    = options.max_compute;
    end
    if isfield(options,'ub') == 1
        ub    = options.ub;
        if length(ub) ~=  length(prior.tau)
            error('Mismatch between the size of upper bounds and the param vector');
        end
    end
    if isfield(options,'lb') == 1
        lb    = options.lb;
        if length(lb) ~=  length(prior.tau)
            error('Mismatch between the size of lower bounds and the param vector');
        end
    end
end

% Confidence Interval Points and Index for confidence serts
alpha  = 1 - conf_sig;
talpha = abs(tinv(alpha/2,T-ny-1));
sort_idx   = round((0.5 + [-conf_sig, conf_sig, 0]/2) * K);

%**************************************************
% Construct the RHS matrix
%**************************************************

X1 = lagX(y,1:lags);
positions_nylags = 1  : size(X1,2);
nylags           = length(positions_nylags);

% forecast launching point
fdata_initval            = lagX(y(end-lags+1:end, :),0:lags-1);
fdata_initval(1:end-1,:) = [];
% add controls if any
if controls_ == 1
    a  = size(X1,2) + 1;
    X1 = [X1 lagX(controls,1:lags)];
    b  = size(X1,2);
    position_controls = a:b;
    % forecast launching point
    lastvalc = lagX(controls(end-lags+1:end, :),0:lags-1);
    fdata_initval  = [fdata_initval, controls(end,:)];
    
else
    position_controls = [];
end
% add proxy shocks if any
if proxy_ == 1
    a     = size(X1,2) + 1;
    X1    = [X1 proxy];
    b     = size(X1,2);
    position_proxy = a:b ;
    % forecast launching point (assume shocks are zero)
    fdata_initval  = [fdata_initval,  zeros(1,ns)];
else
    position_proxy = [];
end
X_ = X1;
if noconstant == 0
    X_ = [X1 ones(T,1)];
    position_constant = size(X_,2);
    % forecast launching point
    fdata_initval  = [fdata_initval, 1];
else
    position_constant = [];
end
nx           = 1 - noconstant;
% other checks
if (lags + hor) >= size(X_,1)
    error('More parameters than observations: consider reducing ''hor'' or ''lags''')
end


%**************************************************
% pre allocation
%**************************************************
Omegaproxy  = eye(ny);
Omega       = eye(ny);
ir_lp       = nan(ny,hor+1,ny,3);  % variable, horizon, shock and mean upper lower
irproxy_lp  = nan(ny,hor+1,1,3);
forecasts   = nan(hor,ny,3);
if dummy == 2
    ir_blp                  = nan(ny,hor+1,ny,K);
    irproxy_blp             = nan(ny,hor+1,1,K);
    bforecasts_no_shocks    = nan(hor,ny,K);
    bforecasts_with_shocks  = nan(hor,ny,K);
    log_dnsty               = nan(hor,1);
end

%**************************************************
%* Computing the LP and DF
%**************************************************
wb = waitbar(0, 'Direct Methods');
for hh = 0 : hor % iteration over horizon
    ytmp = lagX(y, -hh);
    % Reduced Form estimations
    if robust_se_ ~= 0 % Robust SE
        %  options.robust_se_ = robust_se_;
        options.L     = lags + hh + 1;
        olsreg_(hh+1) = ols_reg(ytmp, X_, options);
    else
        olsreg_(hh+1) = ols_reg(ytmp, X_);
    end
    % Proxy IV identifications
    if proxy_
        irproxy_lp(:, hh+1, :, 2) = olsreg_(hh+1).beta(position_proxy, :)'; % mean
        irproxy_lp(:, hh+1, :, 3) = irproxy_lp(:, hh+1, :, 2) + talpha * olsreg_(hh+1).se(position_proxy, :)'; % upper
        irproxy_lp(:, hh+1, :, 1) = irproxy_lp(:, hh+1, :, 2) - talpha * olsreg_(hh+1).se(position_proxy, :)'; % lower
        if hh == 0
            Omegaproxy(:,1) = olsreg_(hh+1).beta(position_proxy, :)';
            Omegasproxy(:,1,:) =  repmat(Omegaproxy(:,1), 1, 1, K) + ...
                olsreg_(hh+1).se(position_proxy, :)'.* randn(ny,1,K);
        end
    end
    % Choleski/Rotated identification
    if hh == 0
        Omega = chol(olsreg_(hh+1).Serror,'Lower') * Q;
        % generate uncertainty on S (flat prior)
        Sbar = olsreg_(hh+1).error'*olsreg_(hh+1).error;
        df   = olsreg_(hh+1).N - olsreg_(hh+1).K;
        [~, Omegas] = generateOmegas(Sbar,df,K,Q);
        Omegasort = sort(Omegas,3);
        
        ir_lp(:, hh+1, :, 2) = Omega; % mean
        ir_lp(:, hh+1, :, 3) = Omegasort(:,:,sort_idx(2)); % UPPER
        ir_lp(:, hh+1, :, 1) = Omegasort(:,:,sort_idx(1)); % LOWER
    else
        [ir] =  iresponse(olsreg_(hh).beta(positions_nylags, :), eye(ny) , 2, Omega);
        ir_lp(:, hh+1, :, 2) = ir(:, 2, :);  % mean
        [ir3] =  iresponse(olsreg_(hh).beta(positions_nylags, :) + talpha * olsreg_(hh).se(positions_nylags, :), eye(ny) , 2, Omega);
        ir_lp(:, hh+1, :, 3) = ir3(:, 2, :); % upper
        [ir1] =  iresponse(olsreg_(hh).beta(positions_nylags, :) - talpha * olsreg_(hh).se(positions_nylags, :), eye(ny) , 2, Omega);
        ir_lp(:, hh+1, :, 1) = ir1(:, 2, :); % lower
    end
    
    % forecast part
    forecasts(hh+1, :, 2) = (fdata_initval * olsreg_(hh+1).beta);
    forecasts(hh+1, :, 3) = forecasts(hh+1, :, 2) + talpha * diag(olsreg_(hh+1).Serror)' ;
    forecasts(hh+1, :, 1) = forecasts(hh+1, :, 2) - talpha * diag(olsreg_(hh+1).Serror)' ;
    
    %======================================================================
    % Baysian DM
    if dummy == 2 % activating Bayesian Direct methods.        
        if hh == 0
            % cholesky IRF on impact
            ir_blp(:, hh+1, :, :)  = Omegas;
            % proxy IRF on impact
            if proxy_
                irproxy_blp(:, hh+1, :, :) = Omegasproxy;
            end
            % one-step ahead forecast
            bforecasts_no_shocks(hh+1, :, :)   = repmat(forecasts(hh+1, :, 2),1,1,K);
            fnoise = Omega * randn(ny,K);
            bforecasts_with_shocks(hh+1, :, :) = bforecasts_no_shocks(hh+1, :, :) + reshape(fnoise,1,ny,K);
            
            % computing the companion matrix
            F       = [prior.Phi.mean(1 : ny * lags, :)'; eye(ny*(lags-1), ny*lags)];
            % constant
            Fo        = [prior.Phi.mean(end, :)'; zeros(ny * (lags-1), 1)];
            % Shocks Companion
            G       = eye(ny * lags, ny);
            
        else % hh > 0
            %********************************************************
            % Conjugate Prior: MN-IW
            %********************************************************
            if max_prior_tau_ == 1 % Maximize the shrinkage on the VAR coefficients 
                try
                    disp(['***********************************************'])
                    disp(['***********************************************'])
                    disp(['Optimization at horizon ' num2str(hh)])
                    x0   = log(prior.tau(hh));
                    switch max_compute
                        case 1
                            optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
                            [xh,fh,~,~,~,~] = ...
                                fminunc('blp_opt_hyperpara',x0,optim_options,...
                                hh,prior,olsreg_(hh),F,G,Fo,positions_nylags,position_constant);
                            %=====================================================================
                        case 2 % constraint
                            % Set default optimization options for fmincon.
                            optim_options = optimset('display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-6);
                            [xh,fh,~,~,~,~,~] = ...
                                fmincon('blp_opt_hyperpara',x0,[],[],[],[],lb(hh),ub(hh),[],optim_options,y,lags,options);
                            %=====================================================================
                        case 3 % Sims
                            crit = 10e-5;
                            nit  = 10e-4;
                            [fh, xh, ~, ~, ~, ~, ~] = ...
                                csminwel('blp_opt_hyperpara',x0,.1*eye(length(x0)),[],crit,nit,hh,prior,olsreg_(hh),F,G,Fo,positions_nylags,position_constant);
                            %=====================================================================
                        case 7 % Matlab's simplex (Optimization toolbox needed).
                            optim_options = optimset('display','iter','MaxFunEvals',30000,'MaxIter',10000,'TolFun',1e-3,'TolX',1e-3);
                            [xh,fh,~,~] = fminsearch('blp_opt_hyperpara',x0,optim_options,y,lags,options);
                    end
                    fprintf('%s = %0.5g\n','Hyper-parameter Mode ',exp(xh))
                    fprintf('%s = %0.5g\n','Marginal Likelihood ',-fh)
                catch
                    warning('Maximization NOT Successful')
                    disp('Using hyper parameter default values')
                    xh = x0;
                end
                prior.tau(hh) = exp(xh);
                disp(['***********************************************'])
            end
            % constructing the posterior moments given the (optimal) shrinkage 
            [posterior_, ~]   = p2p(hh,prior.tau(hh),prior,olsreg_(hh),F,G,Fo,positions_nylags,position_constant);
            % computing the marginal likelihood
            log_dnsty(hh)     = blp_ml(prior.tau(hh),hh,prior,olsreg_(hh),F,G,Fo,positions_nylags,position_constant);
            % Second moments
            S_inv_upper_chol  = chol(inv(posterior_.S));
            XXi_lower_chol    = chol(posterior_.XXi)';
            % number of regressors
            nk = nylags + nx;
            %**************************************************
            %* Generating draws form the Posterior Distribution
            %**************************************************
            for  d =  1 : K %  Gibbs Sampler                
                %======================================================================
                % Inferece: Drawing from the posterior distribution
                % Step 1: draw from the Covariance
                Sigma = rand_inverse_wishart(ny, posterior_.df, S_inv_upper_chol);
                
                % Step 2: given the Covariance Matrix, draw from the AR parameters
                Sigma_lower_chol = chol(Sigma)';
                Phi1 = randn(nk * ny, 1);
                Phi2 = kron(Sigma_lower_chol , XXi_lower_chol) * Phi1;
                Phi3 = reshape(Phi2, nk, ny);
                Phi  = Phi3 + posterior_.PhiHat;
                
                % Step 3: compute IRF
                blp  =  iresponse(Phi, eye(ny) , 2, Omega); 
                ir_blp(:, hh+1, :, d) = blp(:, 2, :);  % mean
                if proxy_
                    blpproxy  =  iresponse(Phi, eye(ny) , 2, Omegaproxy);
                    irproxy_blp(:, hh+1, :, d) = blpproxy(:,2,1);
                end
                % Step 4: compute Forecasts
                bforecasts_no_shocks(hh+1, :, d) = (fdata_initval(1, [positions_nylags position_constant]) * Phi);
                bforecasts_with_shocks(hh+1, :, d) = (fdata_initval(1, [positions_nylags position_constant]) * Phi) + (Sigma_lower_chol * randn(ny,1))';
            end
        end
    end    
    waitbar(hh/hor, wb);
end
close(wb);

%********************************************************
%* Storing the resutls
%*******************************************************
% dm.olsreg_     = olsreg_;
dm.forecasts   = forecasts;
dm.ir_lp       = ir_lp;
dm.irproxy_lp  = irproxy_lp;
if dummy == 2
    dm.ir_blp                 = ir_blp;
    dm.irproxy_blp            = irproxy_blp;
    dm.bforecasts.no_shocks   = bforecasts_no_shocks;         % trajectories of forecasts without shocks
    dm.bforecasts.with_shocks = bforecasts_with_shocks;       % trajectories of forecasts with shocks
    dm.logmlike               = log_dnsty;
    dm.prior                  = prior;
else
    dm.ir_blp                 = [];
    dm.irproxy_blp            = [];
    dm.bforecasts.no_shocks   = [];
    dm.bforecasts.with_shocks = [];
    dm.logmlike               = [];
    dm.prior                  = [];
end



%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
function [S,Omegas] = generateOmegas(Sbar,df,K,Q)
% generate uncertainty on S
ny      = size(Sbar,1);
S       = nan(ny,ny,K);
Omegas  = nan(ny,ny,K);

S_inv_upper_chol = chol(inv(Sbar));

for  d = 1:K
    S(:,:,d)      = rand_inverse_wishart(ny, df, S_inv_upper_chol);
    Omegas(:,:,d) = chol(S(:,:,d),'Lower') * Q;
end
