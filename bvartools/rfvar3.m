function var=rfvar3(ydata,lags,xdata,breaks,lambda,mu,ww,plr)
%RFVAR3 Reduced-form VAR by least squares, with dummy observations for priors.
%
%   var = RFVAR3(ydata,lags,xdata,breaks,lambda,mu)
%   var = RFVAR3(ydata,lags,xdata,breaks,lambda,mu,ww,plr)
%
%   Inputs
%     ydata    T x nvar   dependent variables; the first lags rows are
%                         initial conditions; prior dummies go at the end
%     lags     scalar     number of lags
%     xdata    T x nx     exogenous variables (same length as ydata), or []
%     breaks   column vector of rows after which there is a break. Allows
%              discontinuities in the data and the addition of dummy
%              observations: one dummy observation is lags+1 rows, separated
%              from the rest by a break. The first lags rows after each
%              break are treated as initial conditions.
%     lambda   weight on the co-persistence dummy (sum of all y stable at
%              their initial levels stays there; 5 is a reasonable start).
%              With lambda < 0 the constant is left out of the dummy.
%     mu       weight on the own-persistence dummies (one per variable; 2 is
%              a reasonable start)
%     ww       T x 1   (optional) rescaling of each observation, for
%              heteroskedasticity of known form (Lenza and Primiceri, 2022);
%              rows of y and X are divided by ww
%     plr      (optional) Prior for the Long Run (Giannone, Lenza and
%              Primiceri, 2019): struct with fields H, phi and ybar; the
%              dummy rows are built here. Empty or omitted: no PLR dummies.
%
%   Output
%     var.B    (nvar*lags+nx) x nvar   coefficients, lag-1 block first
%     var.u    residuals, including the dummy rows
%     var.xxi  inv(X'X)
%     var.y, var.X   data matrices including the dummy rows
%
%   If pre-sample x's are not available, repeating xdata(lags+1,:) or
%   copying xdata(lags+1:2*lags,:) into xdata(1:lags,:) are reasonable
%   substitutes; they are used in the persistence priors.
%
%   See also VARPRIOR, PANDEMIC_VARPRIOR, MATRICTINT, BVAR_ML.

%   Original file by Christopher A. Sims, downloaded from
%   http://sims.princeton.edu/yftp/VARtools/matlab/rfvar3.m; the ww and plr arguments are additions.

if nargin<7 || isempty(ww)
   scale_ = 0;
else
    % correct for heteroskedasticity
    scale_=1;
end
if nargin < 8, plr = []; end

[T,nvar] = size(ydata);
nox = isempty(xdata);
if ~nox
    [T2,nx] = size(xdata);
else
    T2 = T;
    nx = 0;
    xdata = zeros(T2,0);
end
% note that x must be same length as y, even though first part of x will not be used.
% This is so that the lags parameter can be changed without reshaping the xdata matrix.
if T2 ~= T, error('Mismatch of x and y data lengths'),end
if nargin < 4
    nbreaks = 0;
    breaks = [];
else
    nbreaks = length(breaks);
end
breaks = [0;breaks;T];
smpl = [];
for nb = 1:nbreaks+1
    smpl = [smpl;[breaks(nb)+lags+1:breaks(nb+1)]'];
end
Tsmpl = size(smpl,1);
X = zeros(Tsmpl,nvar,lags);
for is = 1:length(smpl)
    X(is,:,:) = ydata(smpl(is)-(1:lags),:)';
end
X = [X(:,:) xdata(smpl,:)];
y = ydata(smpl,:);

% rescale if heteroskedasticity corrected
if scale_ ==1
    ww = [ww; ones(length(y)-length(ww),1) ];
    y = y./ repmat(ww,1,size(y,2)) ;
    X = X ./ repmat(ww,1,size(X,2));    
end
% Everything now set up with input data for y=Xb+e 

% prior dummy rows follow the first break (none in OLS calls); lambda changes
% sign below, so keep it for the diffuse prior on the constant
dummies = nbreaks > 0 && breaks(2) < T;
lambda0 = lambda;

% Add persistence dummies
if lambda ~= 0 || mu > 0
    ybar = mean(ydata(1:lags,:),1);
    if ~nox
        xbar = mean(xdata(1:lags,:),1);
    else
        xbar = [];
    end
    if lambda ~= 0
        if lambda>0
            xdum = lambda*[repmat(ybar,1,lags) xbar];
        else
            lambda = -lambda;
            xdum = lambda*[repmat(ybar,1,lags) zeros(size(xbar))];
        end
        ydum = zeros(1,nvar);
        ydum(1,:) = lambda*ybar;
        y = [y;ydum];
        X = [X;xdum];
    end
    if mu>0
        xdum = [repmat(diag(ybar),1,lags) zeros(nvar,nx)]*mu;
        ydum = mu*diag(ybar);
        X = [X;xdum];
        y = [y;ydum];
    end
end

% Diffuse prior on the deterministic terms (GLP 2019, App. B.2): only the
% co-persistence dummy with lambda > 0 involves the constant (the mu and PLR
% rows have a zero constant column); the row is added only to the
% deterministic columns with no prior information (e.g. not to the pandemic
% dummies, which have their own prior)
if dummies && lambda0 <= 0 && nx > 0
    Vc  = 10e6;
    Xp  = X(breaks(2)-lags+1:end, nvar*lags+1:end);   % prior rows, deterministic columns
    nop = find(all(Xp == 0, 1));                      % columns with no prior information
    E   = eye(nx);
    y   = [y; zeros(numel(nop), nvar)];
    X   = [X; zeros(numel(nop), nvar * lags), E(nop,:) / sqrt(Vc)];
end

% Add Priors-for-the-Long-Run dummies (Giannone, Lenza & Primiceri 2019, JASA,
% 114:526, eq. 10): one row per long-run direction (row of H), dependent value
% and all lags equal to w_i = (H(i,:)*ybar'/phi(i))*inv(H)(:,i)', constant and
% exogenous columns 0 -- same mechanics as the persistence dummies above. With
% H = eye(nvar) and phi = 1/mu these are the rows of the mu block. plr.ybar is
% computed once (bvar_.m, or bvar_ml.m when tuning phi), so the prior and the
% posterior calls append identical rows.
if ~isempty(plr)
    H    = plr.H;
    phi  = plr.phi(:);
    ybar = plr.ybar(:)';                   % y_bar_0 (paper eq. 7), 1 x nvar
    % inv(H) on the FULL H first: phi_i = Inf switches direction i off -- the
    % row is DROPPED, not zero-padded (a zero row still shifts the dof inside
    % matrictint through gammaln terms that do not cancel; GLP skip likewise).
    Hinv   = H \ eye(nvar);
    active = find(isfinite(phi))';
    m      = numel(active);
    ydum   = zeros(m, nvar);
    xdum   = zeros(m, nvar * lags + nx);
    for j = 1:m
        i          = active(j);
        wi         = (H(i, :) * ybar') / phi(i) * Hinv(:, i)';  % 1 x nvar
        ydum(j, :) = wi;
        xdum(j, :) = [repmat(wi, 1, lags), zeros(1, nx)];       % p lags = w_i, const 0
    end
    y = [y; ydum];
    X = [X; xdum];
end

% Compute OLS regression and residuals
[vl,d,vr] = svd(X,0);
di = 1./diag(d);
B = (vr.*repmat(di',nvar*lags+nx,1))*vl'*y;
u = y-X*B;
xxi = vr.*repmat(di',nvar*lags+nx,1);
xxi = xxi*xxi';

var.B = B;
var.u = u;
var.xxi = xxi;
var.y   = y;
var.X   = X;
