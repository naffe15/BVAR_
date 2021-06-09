function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
%function [ydum,xdum,breaks]=varprior(nv,nx,lags,mnprior,vprior)
% ydum, xdum:   dummy observation data that implement the prior
% breaks:       vector of points in the dummy data after which new dummy obs's start
%                   Set breaks=T+[0;breaks], ydata=[ydata;ydum], xdum=[xdata;xdum], where 
%                   actual data matrix has T rows, in preparing input for rfvar3
% nv,nx,lags: VAR dimensions
% mnprior.tight:Overall tightness of Minnesota prior
% mnprior.decay:Standard deviations of lags shrink as lag^(-decay)
% vprior.sig:   Vector of prior modes for diagonal elements of r.f. covariance matrix
% vprior.w:     Weight on prior on vcv.  1 corresponds to "one dummy observation" weight
%                   Should be an integer, and will be rounded if not.  vprior.sig is needed
%                   to scale the Minnesota prior, even if the prior on sigma is not used itself.
%                   Set vprior.w=0 to achieve this.
% Note:         The original Minnesota prior treats own lags asymmetrically, and therefore
%                   cannot be implemented entirely with dummy observations.  It is also usually
%                   taken to include the sum-of-coefficients and co-persistence components
%                   that are implemented directly in rfvar3.m.  The diagonal prior on v, combined
%                   with sum-of-coefficients and co-persistence components and with the unit own-first-lag
%                   prior mean generates larger prior variances for own than for cross-effects even in 
%                   this formulation, but here there is no way to shrink toward a set of unconstrained 
%                   univariate AR's.

% Original file downloaded from:
% http://sims.princeton.edu/yftp/VARtools/matlab/varprior.m

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


