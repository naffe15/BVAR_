function [dm] = directmethods(y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'directmethods' generates local projeciton IRF and direct forecasts

% Core Inputs:
% - y, data columns variables
% - lags, lag order of the VAR

% Additonal Inputs collected options:

% Filippo Ferroni, 27/02/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Control random number generator
rng('default');
rng(999);

hor          = 24;
conf_sig     = 0.9;
controls_    = 0;
robust_se_   = 1; % by default robust SE
[T,ny]       = size(y);
proxy_       = 0;
noconstant   = 0;
ns           = 0;
Q            = eye(ny);
dummy        = 0;
K            = 5000;
max_prior_tau_ = 0;

% options
if nargin>2
    %======================================================================
    % Various options
    %======================================================================
    if isfield(options,'hor') ==1
        hor = options.hor;
    end
    if isfield(options,'conf_sig') ==1
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
        %     if   (isfield(options,'prior')==1 && strcmp(options.prior,'Phi')==1) || ...
        %             (isfield(options,'prior')==1 && strcmp(options.prior,'Sigma')==1) || ...
        %             (isfield(options,'prior')==1 && strcmp(options.prior,'lambda')==1)
        
        dummy = 2;
        %  warning('The Conjugate prior is still under construction ... ');
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
                if length(prior.Phi.cov) ~= (ny*lags+(1-noconstant) ) && size(prior.Phi.cov,1 )~=size(prior.Phi.cov,2)
                    error('Size mismatch')
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
        elseif isfield(options.priors,'tau') == 1 && strmatch(options.priors.tau,'max')==1 
            max_prior_tau_ = 1;
            prior.tau = ones(hor,1);
        else
            warning(['You did not provide overall shrinkage. Assume to be one at each horizon'])
            prior.tau = 2*ones(hor,1);
        end
    end
    
end

% Confidence Interval
% t(alpha/2,T-Nnar-1)
alpha  = 1 - conf_sig;
talpha = abs(tinv(alpha/2,T-ny-1));

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

wb = waitbar(0, 'Direct Methods');
for hh = 0 : hor
    ytmp = lagX(y, -hh);
    if robust_se_ ~= 0 % Robust SE
        %  options.robust_se_ = robust_se_;
        options.L     = lags + hh + 1;
        olsreg_(hh+1) = ols_reg(ytmp, X_, options);
    else
        olsreg_(hh+1) = ols_reg(ytmp, X_);
    end
    % Proxy IV identification
    if proxy_
        irproxy_lp(:, hh+1, :, 2) = olsreg_(hh+1).beta(position_proxy, :)'; % mean
        irproxy_lp(:, hh+1, :, 3) = irproxy_lp(:, hh+1, :, 2) + talpha * olsreg_(hh+1).se(position_proxy, :)'; % upper
        irproxy_lp(:, hh+1, :, 1) = irproxy_lp(:, hh+1, :, 2) - talpha * olsreg_(hh+1).se(position_proxy, :)'; % lower
        if hh == 0
            Omegaproxy(:,1) = olsreg_(hh+1).beta(position_proxy, :)';
        end
    end
    % choleski identification
    if hh == 0
        %         tmp1_ = ols_reg(ytmp,XX);
        Omega = chol(olsreg_(hh+1).Serror,'Lower') * Q;
        ir_lp(:, hh+1, :, 2) = Omega; % mean
        ir_lp(:, hh+1, :, 3) = Omega; %ir_lp(:,hh+1,:,2) + talpha * 1;%tempoutput.se(position_shock : end,:); % upper
        ir_lp(:, hh+1, :, 1) = Omega; %ir_lp(:,hh+1,:,2) - talpha * 1;%tempoutput.se(position_shock : end,:); % lower
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
    
    
    if dummy == 2 % activating Bayesian Direct methods.
        
        if hh == 0
            % cholesky IRF on impact
            ir_blp(:, hh+1, :, :)  = repmat(Omega,1,1,K);
            % proxy IRF on impact
            if proxy_
                irproxy_blp(:, hh+1, :, :) = repmat(Omegaproxy(:,1),1,1,K);
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
            % Conjugate Prior: N-IW
            %********************************************************
            % constructing the prior mean
            if max_prior_tau_ == 1
                keyboard;
            else
                [posterior,prior] = p2p(hh,prior,olsreg_(hh),F,G,Fo,positions_nylags,position_constant);
                hyperpara         = prior.tau(hh);
                log_dnsty(hh)     = blp_ml(hyperpara,hh,prior,olsreg_(hh),F,G,Fo,positions_nylags,position_constant);
            end
            S_inv_upper_chol    = chol(inv(posterior.S));
            XXi_lower_chol      = chol(posterior.XXi)';
            
            nk = nylags + nx;
            %  Gibbs Sampler
            for  d =  1 : K
                
                %======================================================================
                % Inferece: Drawing from the posterior distribution
                % Step 1: draw from the Covariance
                Sigma = rand_inverse_wishart(ny, posterior.df, S_inv_upper_chol);
                
                % Step 2: given the Covariance Matrix, draw from the AR parameters
                Sigma_lower_chol = chol(Sigma)';
                Phi1 = randn(nk * ny, 1);
                Phi2 = kron(Sigma_lower_chol , XXi_lower_chol) * Phi1;
                Phi3 = reshape(Phi2, nk, ny);
                Phi  = Phi3 + posterior.PhiHat;
                
                % computing IFR
                blp  =  iresponse(Phi, eye(ny) , 2, Omega);
                ir_blp(:, hh+1, :, d) = blp(:, 2, :);  % mean
                if proxy_
                    blpproxy  =  iresponse(Phi, eye(ny) , 2, Omegaproxy);
                    irproxy_blp(:, hh+1, :, d) = blpproxy(:,2,1);
                end                
                % computing Forecasts
                bforecasts_no_shocks(hh+1, :, d) = (fdata_initval(1, [positions_nylags position_constant]) * Phi);
                bforecasts_with_shocks(hh+1, :, d) = (fdata_initval(1, [positions_nylags position_constant]) * Phi) + (Sigma_lower_chol * randn(ny,1))';
               
            end
        end
    end
    
    waitbar(hh/hor, wb);
    % clear ytmp olsreg_
end
close(wb);

% store
% dm.olsreg_     = olsreg_;
dm.forecasts   = forecasts;
dm.ir_lp       = ir_lp;
dm.irproxy_lp  = irproxy_lp;
if dummy == 2
    dm.ir_blp        = ir_blp;
    dm.irproxy_blp   = irproxy_blp;
    dm.bforecasts.no_shocks   = bforecasts_no_shocks;         % trajectories of forecasts without shocks
    dm.bforecasts.with_shocks = bforecasts_with_shocks;       % trajectories of forecasts with shocks
    dm.log_dnsty     = log_dnsty;
else
    dm.irproxy_blp   = [];
    dm.ir_blp        = [];
    dm.bforecasts    = [];
end