function [logL,outputkf] = kfilternan(Phi,Sigma,y,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'kfilternan' runs the Kalman filter and smoother with missing values
% References: Durbin and Koopman (2003)

% Inputs:
% - Phi, AR parameters of the VAR
% - Sigma, Covariance matrix of the reduced form VAR shocks
% - y, data
% - options, various options
% TVP parameters allowed (3rd dimension of Phi,Sigma and options.tauVec
% controls the time of the variation)

% outputkf : (see below)

% State Space Model
% x(t) = A x(t-1) + B Sigma' u(t) ~ N(0,I)
% y(t) = C*(cons + x(t-1))

% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018
% Revised, 9/11/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% kfilternan computes the forward Kalman filter with NaN

[T,var]         = size(y);
data            = y';
tauVec          = ones(T,1);
initialCond     = 0;
adjustment      = 0;
index           = zeros(var,1);
noprint         = 0;
state_space_model = 1; % default VAR state space model
only_logL       = 0;
start           = 1;

if nargin > 3
    if isfield(options,'tauVec') == 1
        tauVec=optinos.tauVec;
    end
    if isfield(options,'initialCond') ==1
        initialCond = options.initialCond;
    end
    if isfield(options,'initialCond')==1 && options.initialCond==2
        initialCond = options.initialCond;
        pZero   = options.pZero;
        aZero   = options.aZero;
    end
    if isfield(options,'adjustment')==1
        adjustment = options.adjustment;
    end
    if isfield(options,'index')==1
        index = options.index;
    end
    if isfield(options,'noprint')==1
        noprint = options.noprint;
    end
    if isfield(options,'state_space_model')==1
        state_space_model = options.state_space_model;
    end
    if isfield(options,'only_logL')==1 % stop at computing the likelihood
        only_logL = options.only_logL;
    end
    if isfield(options,'start')==1 % likelood computed from start until end, default start =1
         start = options.start;
    end
end

if length(tauVec)~=T
    warning('Tauvec is Longer than T');
    quer('c');
end

%==========================================================================
% 1.0 obtain the steady state representation 
switch state_space_model
    case 1
        %------------------------------------------------------------------
        % VAR
        for jj  = 1 : size(Phi,3)
            [A(:,:,jj),B(:,:,jj),C(:,:,jj),const(:,jj),Sigma(:,:,jj),~,index_var]=var2ss(Phi,Sigma,index);
        end
    case 2
        %------------------------------------------------------------------
        % Unobserved Component I(2)
        [A,B,C,const,Sigma,index_var]=uc2ss(Phi,Sigma);
    case 3
        %------------------------------------------------------------------
        % Dynamic Factor Model 
        % TBA
        % [A,B,C,const,Sigma]=dfm2ss(Phi,Sigma,nfac);
end


%==========================================================================
% 1.1. Dimensions and storage
ns      = size(B,1);
vt      = zeros(var,T);
finvt   = zeros(var,var,T);
kpartg  = zeros(ns,var,T);
logLnc  = zeros(T,1);
yfor    = zeros(size(C,1),T);
sfor    = zeros(size(A,1),T);
% Matrices with one additional entrdataStru.data (initialization)
% to recover observables
stt        = zeros(ns,T+1);
ptt        = zeros(ns,ns,T+1);
W          = eye(var);
Zdim       = zeros(T,1);
mat_obspos = zeros(T,var);

% 1.2 Initialization
if initialCond==0
    stt(:,1) = zeros(ns,1);
    P0 = lyapunov_symm(A(:,:,tauVec(1)),...
        B(:,:,tauVec(1))*(Sigma(:,:,tauVec(1))')...
        *Sigma(:,:,tauVec(1))*(B(:,:,tauVec(1))'));
    ptt(:,:,1)=P0;
elseif initialCond==1 %non stationary data
    P0         = 10*eye(size(A,1));
    stt(:,1)   = zeros(size(A,1),1);
    ptt(:,:,1) = P0;
elseif initialCond==2
    P0         = pZero;
    stt(:,1)   = aZero;
    ptt(:,:,1) = P0;
end

nbreak = 0;
time   = 0;
state  = stt(:,1);

%==========================================================================
% 1.3 Start Forward Filter using KF_DK
for ii=1:T
    
    % Handling of missing observations
    ytt     = data(:,ii);
    
    % Determine W and position of the NAN
    ind         = ~isnan(ytt);
    rowt        = find(~isnan(ytt));
    ytt         = ytt(ind);
    Zdim(ii)    = length(ytt);
    Ztt         = W((ind==1),:)*C(:,:,tauVec(ii));
    mat_obspos(ii,1:Zdim(ii)) = rowt;
    dimt        =( 1:Zdim(ii) );
    % Demeaning is done here
    ytt = ytt - Ztt*const(:,tauVec(ii));
    
    % Forecast Part
    % sfor is the state at time t conditional on info at time t-1, s(t|t-1)
    sfor(:,ii)      = A(:,:,tauVec(ii))*state;
    yfor(:,ii)      = C(:,:,tauVec(ii))*(A(:,:,tauVec(ii))*state + const(:,tauVec(ii)));
    
    %     % computing the 1, 2,3,4 step ahead forecast
    %     yfrsct(ii,dimt,1) = yfor(dimt,ii);
    %     yfrsct(ii,dimt,2) = Ztt*(G(:,:,tauVec(ii))^2*stt(:,ii) + C(:,tauVec(ii)));
    %     yfrsct(ii,dimt,3) = Ztt*(G(:,:,tauVec(ii))^3*stt(:,ii) + C(:,tauVec(ii)));
    %     yfrsct(ii,dimt,4) = Ztt*(G(:,:,tauVec(ii))^4*stt(:,ii) + C(:,tauVec(ii)));
    
    [stt(:,ii+1),ptt(:,:,ii+1),logLnc(ii),vt(dimt,ii),finvt(dimt,dimt,ii),...
        kpartg(:,dimt,ii),]=feval(@kf_dk,ytt,Ztt,...%Z(:,:,tauVec(ii)),...
        state,ptt(:,:,ii),A(:,:,tauVec(ii)),...
        B(:,:,tauVec(ii))*(Sigma(:,:,tauVec(ii))'));
    
    % if there is break
    if ii < T
        if tauVec(ii+1) - tauVec(ii) > 0
            nbreak              = nbreak +1;
            time(nbreak)        = ii+1;
            state = stt(:,ii+1) + adjustment(:,nbreak);
        else
            state = stt(:,ii+1);
        end
    end
    
end

%==========================================================================
% 2. Likelihood with Integration Constant
outputkf.logLncFull     = logLnc;
%logLnc                  = logLnc(dataStru.trainVec(1):dataStru.trainVec(2));
% logLnc                  = logLnc(1:end);
logL                    = -0.5*sum(Zdim (start:end) )*log(2*pi)+sum(logLnc (start:end));
if only_logL == 1
    return;
end

%==========================================================================
% 3. Truncate filters and obtain initial observations
% outSt.yferr=vt';
outputkf.yferr = (data- yfor)';
yfor        = yfor';
stt         = stt(:,2:end);
ptt         = ptt(:,:,2:end);

% % add the 1,2,3,4 ste ahead forecasts in the output
% for hf = 1
%     tmp  = (dataStru.data(:,1+hf:end) - yfrsct(1:end-hf,:,hf) ) *  ...
%         (dataStru.data(:,1+hf:end) - yfrsct(1:end-hf,:,hf) )';
%     outSt.frscterror(:,hf) =
% end

%==========================================================================
% 4. Disturbance smoother with TV matrices
% Obtain the Innovations using a disturbance smoother
if state_space_model ==2
    etamat      = zeros(3,T);
else
    etamat      = zeros(var,T);
end
smooth_st   = zeros(ns,T);
rmat        = zeros(ns,T);

%==========================================================================
% 4.1 Initialize RSTAR & start at t=Nobs
rstar       = zeros(ns,1);
[rstar,etamat(:,end)]   = ...
    smoothdis(rstar,...
    (Sigma(:,:,tauVec(end))')*Sigma(:,:,tauVec(end)),B(:,:,tauVec(end))',...
    Ztt', finvt(dimt,dimt,end),zeros(ns),vt(dimt,end));
smooth_st(:,end)        = stt(:,end)+ptt(:,:,end)*rstar;
rmat(:,end)             = rstar;
Nmat = zeros(ns,ns,T);
%==========================================================================
% 4.2 Begin Backward recursion
for ii=(T-1):-1:1
    
    % Varying dimension in the backward recurions
    dimt=( 1:Zdim(ii) );
    % [Ztt]'= [Wtt*Z]' = = Z'*Wtt'
    Ztt = (C(:,:,tauVec(ii))')*( W( mat_obspos(ii,1:Zdim(ii)),:)');
    
    [rstar,etamat(:,ii)]=smoothdis(rstar,...
        (Sigma(:,:,tauVec(ii))')*Sigma(:,:,tauVec(ii)),...
        B(:,:,tauVec(ii))',Ztt,finvt(dimt,dimt,ii),...
        ((A(:,:,tauVec(ii+1))-A(:,:,tauVec(ii+1))*...
        kpartg(:,dimt,ii)*Ztt')'),vt(dimt,ii));
    
    smooth_st(:,ii)=stt(:,ii)+ptt(:,:,ii)*rstar;
    
    rmat(:,ii)=rstar;
    
end
%==========================================================================
% 4.3 smoothed initial condition
if initialCond~=2
    a0 = P0*A(:,:,tauVec(1))'*rstar;   % Note: this is only correct in the case where a0 = zeros(ns,1)
else
    iGG = pinv(A(:,:,tauVec(1)));
    a0  = iGG*( smooth_st(:,1) - B(:,:,tauVec(1))* etamat(:,1));
end

etamat      = (etamat)';
smooth_st   = (smooth_st)';
stt         = stt';

%==========================================================================
% % 5. Check Smoother
% % Check that Smooth States are identical if using disturbance smoother
% % (above) vs. state smoother (below) and if can also recover the observables
tol=1e-5;

Ydem = data - repmat((C(:,:,tauVec(1))*const(:,tauVec(1))),1,size(data,2));
if max(tauVec) > 1 % only one break
    indx = min(find(tauVec==2));
    tmp1      = data(:,1:indx-1) - repmat((C(:,:,tauVec(1))*const(:,tauVec(1))),1,indx-1);
    tmp2      = data(:,indx:end) - repmat((C(:,:,tauVec(2))*const(:,tauVec(2))),1,size(data,2)-indx+1);
    Ydem = [tmp1 tmp2];
end

maxdifYf     = max(max(abs(stt(:,index_var) - Ydem')));
if noprint,
else
    disp(['Max Discrepancy Filtered vs. Actual Data: ' num2str(maxdifYf)]);
end

maxdifYs     = max(max(abs(smooth_st(:,index_var) - Ydem')));
if noprint,
else
    disp(['Max Discrepancy Smooth vs. Actual Data: ' num2str(maxdifYs)]);
end% maxdifY=comparemat(dataStru.data,Ynanfill');

if maxdifYf > tol || maxdifYs > tol
    warning('Smoother and Filter discrepancy exceeds tolerance')
end


%==========================================================================
% 6. Store results

outputkf.index_var          = index_var;
outputkf.logL               = logL;
outputkf.filteredSt         = stt;
outputkf.CovSt              = ptt;
outputkf.smoothSt           = smooth_st;
outputkf.innovations        = etamat;
outputkf.forecastObs        = yfor;
outputkf.aZero              = a0;
outputkf.adjustment         = adjustment;
outputkf.nbreak             = nbreak;
outputkf.time               = time;
outputkf.smoothSt_plus_ss   = zeros(size(smooth_st));
outputkf.filteredSt_plus_ss = zeros(size(smooth_st));
if time>0
    outputkf.smoothSt_plus_ss(1:time-1,:) = smooth_st(1:time-1,:) + repmat(const(:,1)',time-1,1);
    outputkf.smoothSt_plus_ss(time:end,:) = smooth_st(time:end,:) + repmat(const(:,2)',T-time+1,1);
    outputkf.filteredSt_plus_ss(1:time-1,:) = stt(1:time-1,:) + repmat(const(:,1)',time-1,1);
    outputkf.filteredSt_plus_ss(time:end,:) = stt(time:end,:) + repmat(const(:,2)',T-time+1,1);
    
else
    outputkf.smoothSt_plus_ss = smooth_st + repmat(const(:,end)',size(smooth_st,1),1);
    outputkf.filteredSt_plus_ss = stt + repmat(const(:,end)',size(smooth_st,1),1);
end

if state_space_model ==2
    return;
end

Ydem(find(isnan(Ydem))) = 0;
tmpf = outputkf.yferr';
tmpf(find(isnan(tmpf))) = 0;
outputkf.r2 = 1 - diag(tmpf * tmpf') ./ diag(Ydem*Ydem');
%==========================================================================
% % 7. Simulating the state vector
if state_space_model ==2
    etatilde     = zeros(3,T);
else
    etatilde     = zeros(var,T);
end
for tt =1 : T
    Qt = (Sigma(:,:,tauVec(tt))')*Sigma(:,:,tauVec(tt));
    Ct = Qt - Qt*B(:,:,tauVec(tt))'*Nmat(:,:,tt)*B(:,:,tauVec(tt))*Qt;
    [a,b,~] = svd(Ct);
    iS      = a* sqrt(b);
    if state_space_model ==2
        dt = iS*randn(3,1);
    else
        dt = iS*randn(var,1);
    end
    etatilde(:,tt) = dt + Qt*B(:,:,tauVec(tt))'*rmat(:,tt);
end
if nbreak == 0
    [~,smooth_sim]=kfilterRegSplitSimulation(a0,etatilde);
else
    options.nbreak       = nbreak;
    options.time         = time;
    options.adjustment   = adjustment;
    [~,smooth_sim] = kfilterRegSplitSimulation(a0,etatilde,options);
end
if time>0
    outputkf.smoothSt_sim_plus_ss(1:time-1,:) = smooth_sim(1:time-1,:) + repmat(const(:,1)',time-1,1);
    outputkf.smoothSt_sim_plus_ss(time:end,:) = smooth_sim(time:end,:) + repmat(const(:,2)',T-time+1,1);
else
    outputkf.smoothSt_sim_plus_ss = smooth_sim + repmat(const(:,end)',size(smooth_sim,1),1);
end


%% Subroutine kfilterRegSplitSimulation allows to simulate the model
% Inputs
% sInitial: [ns 1] Initial State
% innovMat: [nx T] matrix of innovations
% By being a sub-routine it has access to all variables defined above
% Be-careful not to use an index (ii,jj) used above or to repeat variable
% names
    function [ySim,sSim]=kfilterRegSplitSimulation(sInitial,innovMat,options)
        if nargin < 3
            nbreak = 0;
            adjustment =0;
            time =0;
        else
            if isfield(options,'time')
                time = options.time;
            else
                error('time of the break is missing')
            end
            if isfield(options,'adjustment')
                adjustment = options.adjustment;
            else
                error('State adjustment is missing')
            end
        end
        
        if ~isequal(size(innovMat),[var T])
            error('Input innovMat must be [nx T]')
        end
        sSim=zeros(ns,T);
        ySim=zeros(size(data));
        sSim(:,1)=A(:,:,tauVec(1))*sInitial + ...
            B(:,:,tauVec(1))*innovMat(:,1);
        ySim(:,1)=C(:,:,tauVec(1))*sSim(:,1); %+C(:,tauVec(1)));
        for kk=2:T;
            if any(kk == time)
                % sSim(:,kk-1) = sSim(:,kk-1) ;
                sSim(:,kk-1) = sSim(:,kk-1) + adjustment(:,find(kk == time));
            end
            sSim(:,kk)=A(:,:,tauVec(kk))*sSim(:,kk-1)+...
                B(:,:,tauVec(kk))*innovMat(:,kk);
            ySim(:,kk)=C(:,:,tauVec(kk))*sSim(:,kk); %+C(:,tauVec(kk)));
        end
        ySim=ySim';
        sSim=sSim';
    end

    function [rzero,what]=smoothdis(rone,Q,Rtr,Ztr,Fnvr,Ltr,v)
        % function [rzero,what]=smoothdis(rone,Q,Rtr,Ztr,Fnvr,Ltr,v)
        %
        % Disturbance smoother recursion for a model with no error in the obs.
        % equation using the formulas in Durbin and Koopman Ch. 4.2
        % rzero=Z'*F^(-1)*v + L'*rone
        % what =Q*R'*rone
        %
        % Inputs
        % Rtr = R'
        % Ztr = Z'
        % Ltr = L'
        % Finvr = F^(-1)
        %
        % Revised, 2/15/2017
        % Revised, 3/21/2018
        
        rzero=Ztr*Fnvr*v + Ltr*rone ;
        what=Q*Rtr*rzero;
    end


%% End of File
end

