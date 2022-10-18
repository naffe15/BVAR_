function [CVAR] = cvar_(y,lags,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'cvar_' estimates the paramters of VAR model with classical methods
% Impulse Response Function CI are computed using bootstrap methods

% Core Inputs:
% - y, data columns variables
% - lags, lag order of the VAR

% Additonal Inputs collected options:
% - options are not mandatory. 
% See the Hitchhiker's guide for more details. 
% https://github.com/naffe15/BVAR_/blob/master/HitchhikerGuide_.pdf

% Output: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('the cvar_ funtion needs at least two inputs: data and number of lags');
end
if lags < 1
    error('lags cannot be zero or negative');
end
% number of observable variables
ny                  = size(y, 2);

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
long_run_irf        = 0;            % when 0, it does not compute long run IRF
irf_1STD            = 1;            % when 1, IRF are computed as 1SD increase. Else, IRF are compued as unitary increase in the shock
cfrcst_yes          = 0;            % no conditional forecast unless defined in options
% non_explosive_      = 0;            % 
heterosked          = 0;
replacement         = 1;           % with replacement when bootstrapping
bootstrap           = 1;


signs_irf           = 0;
narrative_signs_irf = 0;
zeros_signs_irf     = 0;
proxy_irf           = 0;
heterosked_irf      = 0;
hmoments_signs_irf  = 0;
nexogenous          = 0;
exogenous           = [];
cnnctdnss_          = 0;
Ridge_              = 0; 
Lasso_              = 0;
ElasticNet_         = 0;
set_irf             = 0;

if any(any(isnan(y))) %== 1
    error('cvar_ cannot handle missing observation; use bvar_ instead');

end

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
%     if isfield(options,'non_explosive_')==1
%         non_explosive_  = options.non_explosive_;
%     end    
    if isfield(options,'heterosked_weights')==1
        ww  = options.heterosked_weights;
        heterosked = 1;
    end    
    %======================================================================
    % Bootstrap options
    %======================================================================
    if isfield(options,'bootstrap')==1 
        bootstrap = options.bootstrap;
    end
    if isfield(options,'replacement')==1
        replacement = options.replacement;
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
    if isfield(options,'hmoments')==1
        if signs_irf  == 0
            warning('You did not provide any sign restrictions.')
            signs{1} = 'isempty(y(1,1,1))==0';
        end
        if signs_irf  == 1
            signs_irf       = 0;  % disactivating signs
        end
        hmoments_signs_irf = 1;
        hmoments           = options.hmoments;
        [f]                = hmoments2matrix(hmoments,ny);
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

    end    
    %======================================================================
    % Regularization options
    %======================================================================    
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
if (ElasticNet_==1 && Lasso_==1 && Ridge_==1) || ...
        (ElasticNet_==1 && Lasso_==1) || ...
        (Lasso_==1 && Ridge_==1) || ...
        (ElasticNet_==1 && Ridge_==1)    
%    fprintf('Incorrect or missing specification of the number of observations. nobs can be at most %4u\n',size(y,1)-firstobs+1);
    error('You chose more than one regularization estimator, Ridge Lasso or ElasticNet. Please chose only one.')
end

%********************************************************
%* Estimate VAR
%********************************************************

idx = firstobs+presample-lags:firstobs+nobs-1;
nx  = 1;
if noconstant
    nx = 0;
end

% organize data as  yy = XX B + E
[yy,XX] = YXB_(y(idx, :),lags,[nx timetrend]);

% preparing the data for LS estimation
ydata   = y(idx, :);
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

% OLS estimate:
if heterosked == 0
    varols  = rfvar3(ydata, lags, xdata, [T; T], 0, 0);    
else
    varols  = rfvar3(ydata, lags, xdata, [T; T], 0, 0, ww);
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

%**************************************************
%* Generating draws form the Posterior Distribution
%**************************************************

% Preallocation of memory
% Matrices for collecting draws from Posterior Density
% Last dimension corresponds to a specific draw
Phi_boots     = zeros(ny*lags+nx+timetrend + nexogenous, ny, K);   % Autoregressive Parameters
Sigma_boots   = zeros(ny,ny,K);                     % Shocks Covariance
ir_boots      = zeros(ny,hor,ny,K);                 % variable, horizon, shock and draws - Cholesky IRF
irlr_boots    = zeros(ny,hor,ny,K);                 % variable, horizon, shock and draws - Long Run IRF
Qlr_boots     = zeros(ny,ny,K);                     % long run impact matrix
e_boots       = zeros(size(yy,1), ny,K);                  % residuals
yhatfut_no_shocks         = NaN(fhor, ny, K);   % forecasts with shocks
yhatfut_with_shocks       = NaN(fhor, ny, K);   % forecast without the shocks
yhatfut_cfrcst            = NaN(fhor, ny, K);   % forecast conditional on endogenous path
if signs_irf == 1
    irsign_boots = ir_boots;
    Omega_boots    = Sigma_boots;
end
if nexogenous > 0 
    irx_boots = zeros(ny,hor,nexogenous,K);  
    % Ox_boots  = Sigma_boots;
end
if narrative_signs_irf == 1
    irnarrsign_boots = ir_boots;
    Omegan_boots     = Sigma_boots;
end
if hmoments_signs_irf == 1
    irhmomsign_boots = ir_boots;
    Omegam_boots     = Sigma_boots;
end
if zeros_signs_irf == 1
    irzerosign_boots   = ir_boots;
    Omegaz_boots       = Sigma_boots;
end
if proxy_irf == 1
    irproxy_boots = ir_boots;
end
if heterosked_irf == 1
   irheterosked_boots = ir_boots;
   Omegah_boots       = Sigma_boots;
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
forecast_data.initval     = ydata(end-lags+1:end, :);

% # of RHS variables
nk               = ny*lags+nx+timetrend + nexogenous;

waitbar_yes = 0;
if K > 99
    waitbar_yes = 1;
    wb = waitbar(0, 'Bootstrapping');
end

% LS quantities for bootstrapping
Phi_    = varols.B;
e_      = varols.u;
Sigma_  = 1/(nobs-nk) * (e_' * e_);

% use penalized LS if activated in options
if Ridge_ == 1
    Phi_    = varRidge.B;
    e_      = varRidge.u;
    Sigma_  = 1/(nobs-nk) * (e_' * e_);
end
if Lasso_ == 1
    Phi_    = varLasso.B;
    e_      = varLasso.u;
    Sigma_  = 1/(nobs-nk) * (e_' * e_);
end
if ElasticNet_ == 1
    Phi_    = varElasticNet.B;
    e_      = varElasticNet.u;
    Sigma_  = 1/(nobs-nk) * (e_' * e_);
end

% setting of the forecast
deterministic_data.xdata = ones(nobs, nx);
if timetrend
    deterministic_data.xdata = [deterministic_data.xdata [1-lags : T-lags]'];
end
if nexogenous>0
    deterministic_data.xdata = [deterministic_data.xdata exogenous(idx,:)];
end    

% Start bootstrap procedure
for  d =  1 : K
    % Bootstrap residuals
    switch bootstrap
        case 1
            % sampling with replacement 
            % generating K new samples by bootstrapping the VAR reduced form errors
            e_boots(:,:,d) = bootstrap_(e_, 1, replacement);
            % generating K new initial conditions
            initval     =  bootstrap_(ydata(1 : firstobs-1,:), 1, replacement);
        case 2
            % Wild bootstrap based on simple distribution (~Rademacher)
            rr = 1-2*(rand(nobs,1)>0.5);
            %T x n randomly choose the sign of the time T shocks (all)
            e_boots(:,:,d) = (e_ .* (rr*ones(1,ny)));
            initval        = ydata(1 : firstobs-1,:);
    end
    % With bootstrapped errors / intial conditions construct new data
    deterministic_data.initval = initval;    
    [~,ystar0]                 = forecasts(deterministic_data,Phi_,Sigma_,nobs,lags,e_boots(:,:,d),1);        
    ystar                      = [initval; ystar0];
    % Construct left and right ahnd side variables for LS methods
    [yys,XXs] = YXB_(ystar(idx, :),lags,[nx timetrend]);
    % correct for heteroskedasticity of known form (ww) if any
    if heterosked == 1
        wws = [ww; ones(length(yys)-length(ww),1) ];
        yys = yys ./ repmat(wws,1,size(y,2)) ;
        XXs = XXs ./ repmat(wws,1,size(X,2));
    end
    % Compute LS estimator
    [vl,d_,vr] = svd(XXs,0);
    di         = 1./diag(d_);
    Phi        = (vr .* repmat(di', nk, 1)) * vl' * yys;
    % Penalized the LS estimator if activated in options
    if Ridge_ == 1
        xxi        = vr .* repmat(di',nk,1);
        iXX        = xxi * xxi';
        Phi = (eye(nk) + Ridge_lambda*iXX)\Phi;
    end
    if Lasso_ == 1
        LassoPhi = nan(size(Phi));
        for vv = 1 : ny
            LassoPhi(:,vv) = lasso(XXs,yys(:,vv),'Lambda',Lasso_lambda);
        end
        Phi = LassoPhi;
    end
    if ElasticNet_ == 1
        ElasticNetPhi = nan(size(Phi));
        for vv = 1 : ny
            ElasticNetPhi(:,vv) = ...
                lasso(XXs,yys(:,vv),'Lambda',ElasticNet_lambda,'Alpha',ElasticNet_alpha);
        end
        Phi = ElasticNetPhi;            
    end
    % compute residual and covariance matrix
    errors  = yys - XXs*Phi;
    Sigma   = 1/(nobs-nk)*(errors' * errors);        
    % store the draw
    Phi_boots(:,:,d)   = Phi;
    Sigma_boots(:,:,d) = Sigma; 
        
    %======================================================================
    % IRF
    % Compute the impulse response functions
    % with cholesky
    if irf_1STD == 1
        % one STD increase
        ir_boots(:,:,:,d)      = iresponse(Phi(1 : ny*lags, 1 : ny),Sigma,hor,eye(ny));
    else
        % one percent increase
        ir_boots(:,:,:,d)      = iresponse(Phi(1 : ny*lags, 1 : ny),Sigma,hor,eye(ny),0);
    end
    if nexogenous > 0
        Phi1         = Phi(1 : ny*lags + nx + timetrend, 1 : ny);
        Exo1         = zeros(ny,nexogenous);
        Exo1(:,1:nexogenous) = Phi(ny*lags + nx + timetrend +1 ....
                            : ny*lags + nx + timetrend + nexogenous, 1 : ny)';
        Sig1         = eye(ny);
        % unitary increase        
        irx_boots(:,:,:,d) = iresponse(Phi1,Sig1,hor,Exo1);        
    end
    % define the identity rotation    
    Omega = eye(ny);
    
    % with long run restrictions
    if long_run_irf == 1
        [irlr,Qlr]             = iresponse_longrun(Phi(1 : ny*lags, 1 : ny),Sigma,hor,lags);
        irlr_boots(:,:,:,d)    = irlr;
        Qlr_boots(:,:,d)       = Qlr;
        Omega                  = Qlr;
    end
%     % with sign restrictions
%     if signs_irf == 1
%         [irsign,Omega]         = iresponse_sign(Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs);
%         irsign_boots(:,:,:,d)  = irsign;
%         Omega_boots(:,:,d)     = Omega;
%     end
%     % with narrative and sign restrictions
%     if narrative_signs_irf == 1
%         [irnarrsign,Omega]         = iresponse_sign_narrative(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,narrative);
%         irnarrsign_boots(:,:,:,d)  = irnarrsign;
%         Omegan_boots(:,:,d)        = Omega;
%     end
%     % with higher-moments and sign restrictions
%     if hmoments_signs_irf == 1
%         [irhmomsign,Omega]         = iresponse_sign_hmoments(errors,Phi(1 : ny*lags, 1 : ny),Sigma,hor,signs,hmoments,f);
%         irhmomsign_boots(:,:,:,d)  = irhmomsign;
%         Omegam_boots(:,:,d)        = Omega;
%     end
    % with zeros and sign restrictions
    if zeros_signs_irf == 1         %= iresponse_zeros_signs( Phi,Sigma,bvar1.hor,lags,var_pos,f,sr);
        [irzerosign,Omega]          = iresponse_zeros_signs(Phi,Sigma,hor,lags,var_pos,f,sr);
        irzerosign_boots(:,:,:,d)   = irzerosign;
        Omegaz_boots(:,:,d)         = Omega;
    end
    % with proxy
    if proxy_irf == 1
        in.res                  = e_boots(:,:,d);
        in.Phi                  = Phi_boots(:,:,d)  ;
        in.Sigma                = Sigma;
        tmp_                    = iresponse_proxy(in);
        irproxy_boots(:,:,1,d)  = tmp_.irs';
        clear tmp_
    end
    % with heteroskedasticity 
    if heterosked_irf == 1
        [irheterosked,Omegah]       = iresponse_heterosked(Phi(1 : ny*lags, 1 : ny),errors,hor,heterosked_regimes);
        irheterosked_boots(:,:,:,d) = irheterosked;
        Omegah_boots(:,:,d)         = Omegah;
    end
    
    %======================================================================
    % Forecasts
    % compute the out of sample forecast (unconditional)
    [frcst_no_shock,frcsts_with_shocks] = forecasts(forecast_data,Phi,Sigma,fhor,lags);
    yhatfut_no_shocks(:,:,d)            = frcst_no_shock;
    yhatfut_with_shocks(:,:,d)          = frcsts_with_shocks;
    
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
    % Connectedness
    if cnnctdnss_ > 0 
%         if use_omega == 1
%             [C] = connectedness(Phi,Sigma,nethor,Omega); 
%         else
%         end
        if cnnctdnss_ == 2 
            Sigma_lower_chol = chol(Sigma)';
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
%     dd = 0; % reset 
end
if waitbar_yes, close(wb); end


%********************************************************
%* Storing the resutls
%*******************************************************

%==========================================================================
% classical inference: LS estimator
CVAR.Phi_ols    = varols.B;
CVAR.e_ols      = varols.u;
CVAR.Sigma_ols  = 1/(nobs-nk)*varols.u'*varols.u;
[CVAR.InfoCrit.AIC, CVAR.InfoCrit.HQIC, CVAR.InfoCrit.BIC] = IC(CVAR.Sigma_ols, CVAR.e_ols, nobs, nk);
% the model with the lowest IC is preferred

% OLS irf
% with cholesky
CVAR.ir_ols      = iresponse(CVAR.Phi_ols(1 : ny*lags, 1 : ny),CVAR.Sigma_ols,hor,eye(ny));
% with long run
if long_run_irf == 1
    [irlr,Qlr]              = iresponse_longrun(CVAR.Phi_ols(1 : ny*lags, 1 : ny),CVAR.Sigma_ols,hor,lags);
    CVAR.irlr_ols           = irlr;
    CVAR.Qlr_ols(:,:)       = Qlr;
end
% set identified IRF
if set_irf > 0  
    wb = waitbar(0, ['Generating rotations for set-identification - LS Estimator']);    
    % with sign restrictions
    if signs_irf == 1
        for d1 = 1 : set_irf
            [CVAR.irsign_ols(:,:,:,d1),CVAR.Omega_ols(:,:,d1)] = ...
                iresponse_sign(CVAR.Phi_ols(1 : ny*lags, 1 : ny),CVAR.Sigma_ols,hor,signs);
            waitbar(d1/set_irf, wb);
        end
    end
    % with narrative and sign restrictions
    if narrative_signs_irf == 1
        for d1 = 1 : set_irf
            [CVAR.irnarrsign_ols(:,:,:,d1),CVAR.Omegan_ols(:,:,d1)] = ...
                iresponse_sign_narrative(CVAR.e_ols,CVAR.Phi_ols(1 : ny*lags, 1 : ny),CVAR.Sigma_ols,hor,signs,narrative);
            waitbar(d1/set_irf, wb);
        end
    end
    % with zeros and sign restrictions
    if zeros_signs_irf == 1         
        for d1 = 1 : set_irf
            [CVAR.irzerosign_ols(:,:,:,d1),CVAR.Omegaz_ols(:,:,d1)] = ...
                iresponse_zeros_signs(CVAR.Phi_ols,CVAR.Sigma_ols,hor,lags,var_pos,f,sr);
            waitbar(d1/set_irf, wb);
        end
    end
    close(wb)
end
% proxy
if proxy_irf == 1
    inols.res               = CVAR.e_ols;
    inols.Phi               = CVAR.Phi_ols;
    inols.Sigma             = CVAR.Sigma_ols;
    inols.compute_F_stat    = 1;
    tmp_                    = iresponse_proxy(inols);
    CVAR.irproxy_ols(:,:,1) = tmp_.irs';
    CVAR.proxy.F_m          = tmp_.F_m;
    CVAR.proxy.F_m_rob      = tmp_.F_m_rob;
    CVAR.proxy.R2adj_m      = tmp_.R2adj_m;
    CVAR.proxy.data         = options.proxy;
    clear tmp_
end
% with heteroskedasticity
if heterosked_irf == 1
    [CVAR.irheterosked_ols,CVAR.Omegah_ols] = ...
        iresponse_heterosked(CVAR.Phi_ols(1 : ny*lags, 1 : ny),CVAR.e_ols,hor,heterosked_regimes);    
end
% test the normality of the ols VAR residuals (matlab stat toolbox needed)
if  exist('kstest') ==2
    for gg = 1 : ny
        [H,Pv] = kstest(CVAR.e_ols(:,gg)/sqrt(CVAR.Sigma_ols(gg,gg)));
        CVAR.HP(gg,:) = [H,Pv];
        %      H = 0 => Do not reject the null hypothesis at the 5% significance
        %      level. 
    end
else
    CVAR.HP = [];
end

%==========================================================================
% Penalized Approaches (Regularization)
penalizationstrn = {'Ridge','Lasso','ElasticNet'};
% loop across penalization approaches
for pp = 1 : 3    
    Phi_ = []; Sigma_ = []; u_ = []; Omega_(:,:,1) = eye(ny);  InfoCrit_ = [];   
    if eval([ penalizationstrn{pp} '_ == 1'])
        % AR param
        eval(['CVAR.' penalizationstrn{pp} '.Phi    = var' penalizationstrn{pp} '.B;']);
        eval(['Phi_   = CVAR.' penalizationstrn{pp} '.Phi;'])       
        % Error term
        eval(['CVAR.' penalizationstrn{pp} '.e      = var' penalizationstrn{pp} '.u;']);
        eval(['u_ = CVAR.' penalizationstrn{pp} '.e;'])
        % Covariance Matrix        
        eval(['CVAR.' penalizationstrn{pp} '.Sigma  = 1/(nobs-nk) * var' penalizationstrn{pp} '.u'' * var' penalizationstrn{pp} '.u;'])        
        eval(['Sigma_ = CVAR.' penalizationstrn{pp} '.Sigma;'])
        % Recursive IRFs
        eval(['CVAR.' penalizationstrn{pp} '.ir      = iresponse(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,eye(ny));']);                
        % info crit
        [InfoCrit_.AIC, InfoCrit_.HQIC, InfoCrit_.BIC] ...
            = IC(Sigma_, u_, nobs, nk);
        eval(['CVAR.' penalizationstrn{pp} '.InfoCrit      = InfoCrit_;']);                                     
        % IRFs with different identification schemes
        % point identification: LR
        if long_run_irf == 1
            eval(['[CVAR.' penalizationstrn{pp} '.irlr,Omega_] = iresponse_longrun(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,lags);'])
        end
        % point identification: with heteroskedasticity
        if heterosked_irf == 1
            eval(['[CVAR.' penalizationstrn{pp} '.irheterosked,Omega_] = iresponse_heterosked(Phi_(1 : ny*lags, 1 : ny),u_,hor,heterosked_regimes);'])
        end
        % set identification: signs
        if set_irf > 0  
            wb = waitbar(0, ['Generating rotations for set-identification - ' penalizationstrn{pp} ' Estimator']);
        end        
        if signs_irf == 1
            for d1 = 1 : set_irf                
                eval(['[CVAR.' penalizationstrn{pp} '.irsign_boots(:,:,:,' num2str(d1) ...
                    '),Omega_(:,:,' num2str(d1) ')] = iresponse_sign(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,signs);'])
                waitbar(d1/set_irf, wb);
            end            
        end
        % set identification: narrative and sign restrictions
        if narrative_signs_irf == 1
            for d1 = 1 : set_irf
                eval(['[CVAR.' penalizationstrn{pp} '.irnarrsign_boots(:,:,:,' num2str(d1) ...
                    '),Omega_(:,:,' num2str(d1) ')] = iresponse_sign_narrative(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,signs,narrative);'])
                waitbar(d1/set_irf, wb);
            end
        end
        % set identification: zeros and sign restrictions
        if zeros_signs_irf == 1
            for d1 = 1 : set_irf
                eval(['[CVAR.' penalizationstrn{pp} 'irzerosign_boots.(:,:,:,' num2str(d1) ...
                    '),Omega_(:,:,' num2str(d1) ')] = iresponse_zeros_signs(Phi_(1 : ny*lags, 1 : ny),Sigma_,hor,lags,var_pos,f,sr);'])
                waitbar(d1/set_irf, wb);
            end
        end
        if set_irf>0, close(wb); end
    end
    if cnnctdnss_ == 1 % default identification (Pesaran and Shin)
        eval(['CVAR.' penalizationstrn{pp} '.Connectedness = connectedness(Phi_(1 : ny*lags, 1 : ny),Sigma_,nethor);']);        
    elseif cnnctdnss_ == 2 % customized identification
        Sigma_lower_chol = chol(Sigma_)';
        % this part is not correct as the median orthogonal rotation is not
        % orthogonal
        Omegam           = median(Omega_,3);
        eval(['CVAR.' penalizationstrn{pp} '.Connectedness = connectedness(Phi_(1 : ny*lags, 1 : ny), Sigma_, nethor, Sigma_lower_chol * Omegam);']);                
    end

end
    
%==========================================================================
% Store the output
% the last dimension of these objects corresponds to a bootstrapped sample

% inference and IRFs
CVAR.Phi_boots    = Phi_boots;          % draws from the autoregressive part
CVAR.Sigma_boots  = Sigma_boots;        % draws from the covarance matrix
CVAR.ir_boots     = ir_boots;           % draws from the IRF with cholesky
CVAR.irlr_boots   = irlr_boots;         % draws from the IRF with Long Run
CVAR.Qlr_boots    = Qlr_boots;          % Long Run Rotation matrix
CVAR.lags         = lags;               % lags
CVAR.N            = ny;                 % number of variables
CVAR.e_boots      = e_boots;            % residuals
CVAR.XX           = XX;                 % regressors (no dummy)
CVAR.yy           = yy;                 % dependent  (no dummy)

% prediction
CVAR.fhor         = fhor;               % forecast horizon
CVAR.hor          = hor;                % IRF horizon
CVAR.forecasts.no_shocks      = yhatfut_no_shocks;         % trajectories of forecasts without shocks
CVAR.forecasts.with_shocks    = yhatfut_with_shocks;       % trajectories of forecasts with shocks
CVAR.forecasts.conditional    = [];            % trajectories of conditional forecasts
CVAR.forecasts.EPScond        = [];            % shocks of conditional forecasts
if cfrcst_yes ~= 0
    CVAR.forecasts.conditional    = yhatfut_cfrcst;            % trajectories of forecasts
    CVAR.forecasts.EPScond        = EPS;                       % shocks of forecasts
end
CVAR.forecast_data            = forecast_data;

%
CVAR.varnames     = varnames;
CVAR.ndraws       = K;

if signs_irf == 1 && narrative_signs_irf == 0
    CVAR.irsign_boots = irsign_boots;
    CVAR.Omegas       = Omega_boots;
else
    CVAR.irsign_boots = [];
    CVAR.Omegas       = [];
end
if narrative_signs_irf == 1
    CVAR.irnarrsign_boots = irnarrsign_boots;
    CVAR.Omegan           = Omegan_boots;
else
    CVAR.irnarrsign_boots = [];
    CVAR.Omegan           = [];
end
if hmoments_signs_irf == 1
    CVAR.irhmomsign_boots = irhmomsign_boots;
    CVAR.Omegam           = Omegam_boots;
else
    CVAR.irhmomsign_boots = [];
    CVAR.Omegam           = [];
end
if zeros_signs_irf == 1
    CVAR.irzerosign_boots   = irzerosign_boots;
    CVAR.Omegaz             = Omegaz_boots;
else
    CVAR.irzerosign_boots = [];
    CVAR.Omegaz           = [];
end
if heterosked_irf == 1
    CVAR.irheterosked_boots   = irheterosked_boots;
    CVAR.Omegah_boots         = Omegah_boots;
    CVAR.Omegah               = CVAR.Omegah_boots; 
else
    CVAR.irheterosked_boots   = [];
    CVAR.Omegah_boots         = [];
    CVAR.Omegah               = [];
end
if proxy_irf == 1
    CVAR.irproxy_boots = irproxy_boots;
else
    CVAR.irproxy_boots= [];
end
if nexogenous > 0
    CVAR.irx_boots = irx_boots;
else
    CVAR.irx_boots = [];
end

% connetedness computed using LS (or penalized LS if activated in options)
if cnnctdnss_
    CVAR.Connectedness.Index         = CnndtnssIndex;
    CVAR.Connectedness.FromAlltoUnit = CnndtnssFromAlltoUnit;
    CVAR.Connectedness.FromUnitToAll = CnndtnssFromUnitToAll;
    CVAR.Connectedness.theta         = Ctheta;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end of cvar_.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

