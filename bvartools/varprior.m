function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
%VARPRIOR Dummy observations that implement the Minnesota prior.
%
%   [ydum,xdum,breaks] = VARPRIOR(nv,nx,lags,mnprior,vprior)
%
%   Inputs
%     nv, nx, lags     number of variables, exogenous variables and lags
%     mnprior.tight    overall tightness of the Minnesota prior
%     mnprior.decay    standard deviations of lags shrink as lag^(-decay)
%     vprior.sig       nv x 1   prior modes of the diagonal of the residual
%                               covariance (also used to scale the
%                               Minnesota prior)
%     vprior.w         weight on the covariance prior; 1 is one dummy
%                      observation (rounded to an integer); 0 uses sig only
%                      to scale the Minnesota prior
%
%   Outputs
%     ydum, xdum   dummy observations
%     breaks       rows of the dummy data after which a new dummy
%                  observation starts. To prepare RFVAR3 input with T rows
%                  of data: breaks = T+[0;breaks], ydata = [ydata;ydum],
%                  xdata = [xdata;xdum].
%
%   The original Minnesota prior treats own lags asymmetrically and cannot
%   be implemented entirely with dummy observations. The sum-of-coefficients
%   and co-persistence components are implemented in RFVAR3.
%
%   See also RFVAR3, PANDEMIC_VARPRIOR, BVAR_ML.

%   Original file by Christopher A. Sims, downloaded from
%   http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m

if ~isempty(mnprior)
    xdum = zeros(lags+1,nx,lags,nv);
    ydum = zeros(lags+1,nv,lags,nv);
    for il = 1:lags
        ydum(il+1,:,il,:) = il^mnprior.decay*diag(vprior.sig);
    end
    ydum(1,:,1,:) = diag(vprior.sig);
    ydum = mnprior.tight*reshape(ydum,[lags+1,nv,lags*nv]);
    ydum = flipdim(ydum,1);
    xdum = mnprior.tight*reshape(xdum,[lags+1,nx,lags*nv]);
    xdum = flipdim(xdum,1);
    breaks = (lags+1)*[1:(nv*lags)]';
    lbreak = breaks(end);
else
    ydum = [];
    xdum = [];
    breaks = [];
    lbreak = 0;
end
% elle = 0;
% for j1 = 1 : nv
%     for j2  =  1 : lags
%         elle = 1 + elle;
%         ydum(:,:,elle) = ydum(:,:,elle) * mnprior.unit_root_(j1);
%     end
% end
if ~isempty(vprior) && vprior.w>0
    ydum2 = zeros(lags+1,nv,nv);
    xdum2 = zeros(lags+1,nx,nv);
    ydum2(end,:,:) = diag(vprior.sig);
    for i = 1:vprior.w
        ydum = cat(3,ydum,ydum2);
        xdum = cat(3,xdum,xdum2);
        breaks = [breaks;(lags+1)*[1:nv]'+lbreak];
        lbreak = breaks(end);
    end
end
dimy = size(ydum);
ydum = reshape(permute(ydum,[1 3 2]),dimy(1)*dimy(3),nv);
xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx);
breaks = breaks(1:(end-1));


