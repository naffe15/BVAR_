function var=rfvar3(ydata,lags,xdata,breaks,lambda,mu,ww)
%function var=rfvar3(ydata,lags,xdata,breaks,lambda,mu)
% This algorithm goes for accuracy without worrying about memory requirements.
% ydata:   dependent variable data matrix
% xdata:   exogenous variable data matrix
% lags:    number of lags
% breaks:  rows in ydata and xdata after which there is a break.  This allows for
%          discontinuities in the data (e.g. war years) and for the possibility of
%          adding dummy observations to implement a prior.  This must be a column vector.
%          Note that a single dummy observation becomes lags+1 rows of the data matrix,
%          with a break separating it from the rest of the data.  The function treats the 
%          first lags observations at the top and after each "break" in ydata and xdata as
%          initial conditions. 
% lambda:  weight on "co-persistence" prior dummy observations.  This expresses
%          belief that when data on *all* y's are stable at their initial levels, they will
%          tend to persist at that level.  lambda=5 is a reasonable first try.  With lambda<0,
%          constant term is not included in the dummy observation, so that stationary models
%          with means equal to initial ybar do not fit the prior mean.  With lambda>0, the prior
%          implies that large constants are unlikely if unit roots are present.
% mu:      weight on "own persistence" prior dummy observation.  Expresses belief
%          that when y_i has been stable at its initial level, it will tend to persist
%          at that level, regardless of the values of other variables.  There is
%          one of these for each variable.  A reasonable first guess is mu=2.
%      The program assumes that the first lags rows of ydata and xdata are real data, not dummies.
%      Dummy observations should go at the end, if any.  If pre-sample x's are not available,
%      repeating the initial xdata(lags+1,:) row or copying xdata(lags+1:2*lags,:) into 
%      xdata(1:lags,:) are reasonable subsititutes.  These values are used in forming the
%      persistence priors.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/rfvar3.m

if nargin<7
   scale_ = 0;
else
    % correct for heteroskedasticity
    scale_=1;
end

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
