function log_dnsty = bvar_ml(hyperpara,y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'bvar_ml' computes the marginal likelihood for the minnesota prior 

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

%******************************************************** 
% SETTINGS
%******************************************************** 


% randn('state',999);
% rand('state',999);
% 
ny                  = size(y, 2);
first_obs           = lags+1; 
presample           = 0;
noconstant          = 0;

bvar_prior_tau      = hyperpara(1);%3;
bvar_prior_decay    = hyperpara(2);%0.5;
bvar_prior_lambda   = hyperpara(3);%5;                                    
bvar_prior_mu       = hyperpara(4);%2;
bvar_prior_omega    = hyperpara(5);%1;
unit_root_          = ones(ny,1);
flat                = 0; 
heterosked          = 0;

esse = ones(size(y,1),ny);

if nargin > 3 
%     fields = fieldnames(options);
%     nf     = size(fieldnames(options),1);        
%     for i = 1 : nf
%        if strcmp(fields{i},'bvar_prior_train')
%            bvar_prior_train = options.bvar_prior_train;
%        elseif strcmp(fields{i},'first_obs')
%            first_obs = options.first_obs;  
%        elseif strcmp(fields{i},'presample')
%            presample = options.presample;          
%        elseif strcmp(fields{i},'train')
%            train = options.train;  
%        elseif strcmp(fields{i},'noconstant')
%            noconstant = options.nocostant;  
%        end    
%     end
    if isfield(options,'first_obs')==1
        first_obs = options.first_obs;
    end
    if isfield(options,'presample')==1
        presample = options.presample;
    end
    if isfield(options,'noconstant')==1
        noconstant = options.nocostant;
        %  MINNESOTA PRIOR
    end    
    if isfield(options,'unit_root_')==1
        unit_root_ = options.unit_root_;
        %  MINNESOTA PRIOR: unit root assumption only for a subset of var
    end
    if isfield(options,'esse')==1
        % covid-19 reweighting (see bvar_opt_covid)
        esse = options.esse;
    end
    if isfield(options,'heterosked_weights')==1
        ww  = options.heterosked_weights;
        heterosked = 1;
    end    


end

nobs = size(y,1)-first_obs+1;

if (first_obs+ nobs-1)> size(y,1)
    fprintf('Incorrect or missing specification of the number of observations. nobs can be at most %4u\n',size(y,1)-first_obs+1);
    error('Inconsistent number of observations.') 
end

% Parameters for prior
if first_obs + presample <= lags
    error('first_obs+presample should be > lags (for initializing the VAR)')
end

if first_obs + presample  <= lags
    error('first_obs+presample should be > nlags (for initializating the VAR)')
end

idx = first_obs+presample-lags:first_obs+nobs-1;


if noconstant
        nx = 0;
else
        nx = 1;
end


mnprior.tight = bvar_prior_tau;
mnprior.decay = bvar_prior_decay;
mnprior.unit_root_ = unit_root_;
% Use only initializations lags for the variance prior
vprior.sig            = std(y(first_obs-lags : first_obs+presample-1,:))';
% vprior.sig = std(y(first_obs+presample-lags : first_obs+presample,:))';
vprior.w = bvar_prior_omega;
lambda = bvar_prior_lambda;
mu     = bvar_prior_mu;
[ydum, xdum, pbreaks] = varprior(ny, nx, lags, mnprior, vprior);

ydata = y(idx, :);
T     =  size(ydata, 1);
xdata = ones(T,nx);
% posterior density
if heterosked == 0
    var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu);
else
    var = rfvar3([ydata; ydum], lags, [xdata; xdum], [T; T+pbreaks], lambda, mu, ww);
end
Tu = size(var.u, 1);

% Prior density
Tp = presample + lags;
if nx
    xdata = xdata(1:Tp, :);
else
    xdata = [];
end
varp = rfvar3([ydata(1:Tp, :); ydum], lags, [xdata; xdum], [Tp; Tp + pbreaks], lambda, mu);
Tup = size(varp.u, 1);

prior.df        = Tup - ny*lags - nx - flat*(ny+1);
prior.S         = varp.u' * varp.u;
prior.XXi       = varp.xxi;
prior.PhiHat    = varp.B;

priors.YYdum = varp.y;
priors.XXdum = varp.X;

if prior.df < ny
    error('Too few degrees of freedom in the Inverse-Wishart part of prior distribution. You should increase training sample size.')
end
posterior.df     = Tu - ny*lags - nx - flat*(ny+1);
posterior.S      = var.u' * var.u;
posterior.XXi    = var.xxi;
posterior.PhiHat = var.B; 


%*********************************************************
%* Compute the log marginal data density for the VAR model 
%*********************************************************

posterior_int   = matrictint(posterior.S, posterior.df, posterior.XXi);
prior_int       = matrictint(prior.S, prior.df, prior.XXi);
lik_nobs        = posterior.df - prior.df;    
log_dnsty       = posterior_int - prior_int - 0.5*ny*lik_nobs*log(2*pi);
% minus_log_dnsty = -log_dnsty;
end