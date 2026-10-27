function [ydum,xdum,breaks]=pandemic_varprior(nv,nx,lags,mnprior,vprior, pandemic)
%PANDEMIC_VARPRIOR Minnesota dummy observations plus the Pandemic Priors.
%
%   [ydum,xdum,breaks] = PANDEMIC_VARPRIOR(nv,nx,lags,mnprior,vprior,pandemic)
%
%   Same as VARPRIOR, with xdum widened by pandemic.h columns (one dummy
%   per pandemic period) and a block of dummy observations that shrinks
%   their coefficients. This implements the Pandemic Priors of
%   Cascaldi-Garcia (2025) on top of the Sims-Zha Minnesota prior, rather
%   than on the Banbura, Giannone and Reichlin (2010) version used in the
%   original paper.
%
%   Inputs
%     nv, nx, lags, mnprior, vprior   as in VARPRIOR
%     pandemic.h     number of pandemic periods
%     pandemic.phi   1 x h   shrinkage of the pandemic time dummies
%
%   Outputs
%     ydum, xdum, breaks   as in VARPRIOR (xdum has nx + h columns)
%
%   See also VARPRIOR, RFVAR3, BVAR_.
%
%   Reference: Cascaldi-Garcia, D. (2025), "Pandemic Priors".
%   Built on varprior.m by Christopher A. Sims.

if ~isempty(mnprior)
    xdum = zeros(lags+1,nx + pandemic.h,lags,nv);
    ydum = zeros(lags+1,nv,lags,nv);
    for il = 1:lags
        ydum(il+1,:,il,:) = il^mnprior.decay*diag(vprior.sig);
    end
    ydum(1,:,1,:) = diag(vprior.sig);
    ydum = mnprior.tight*reshape(ydum,[lags+1,nv,lags*nv]);
    ydum = flipdim(ydum,1);
    xdum = mnprior.tight*reshape(xdum,[lags+1,nx + pandemic.h ,lags*nv]);
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
    xdum2 = zeros(lags+1,nx + pandemic.h,nv);
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
xdum = reshape(permute(xdum,[1 3 2]),dimy(1)*dimy(3),nx + pandemic.h);



ydum   = [ydum;  zeros((lags+1)*pandemic.h, nv)];

for j = 1 : pandemic.h
    xnewp  = zeros(lags+1, nx + pandemic.h);
    xnewp(end, nx + j) = pandemic.phi(j);
    xdum   = [xdum;  xnewp];
    pandemic_breaks = size(xdum,1) ;
    breaks          = [breaks; pandemic_breaks];
end


breaks = breaks(1:(end-1));
