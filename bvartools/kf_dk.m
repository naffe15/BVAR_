function   [shatnew,signew,lh,yhat,fin,kgpart,yforc]=kf_dk(y,H,shat,sig,G,M)
%KF_DK One step of the Kalman filter, robust to singular covariances.
%
%   [shatnew,signew,lh,yhat,fin,kgpart,yforc] = KF_DK(y,H,shat,sig,G,M)
%
%   State space:  s(t) = G s(t-1) + M e(t),  e(t) ~ N(0,I);  y(t) = H s(t)
%   (M = R*chol(Q)' when s(t) = G s(t-1) + R n(t) with V(n(t)) = Q).
%
%   Inputs
%     y      ny x 1   observations at t
%     H      ny x ns  observation matrix
%     shat   ns x 1   state estimate at t-1
%     sig    ns x ns  its covariance
%     G      ns x ns  transition matrix
%     M      ns x nm  loading of the shocks
%
%   Outputs
%     shatnew, signew   updated state estimate and covariance
%     lh                log likelihood contribution
%     yhat              forecast error, y - yforc
%     fin               inverse of the forecast error covariance
%                       (generalised inverse)
%     kgpart            partial Kalman gain; the gain is G*kgpart
%     yforc             one-step forecast of y
%
%   See also KFILTERNAN.
%
%   Based on the Kalman filter of Christopher A. Sims, with added outputs
%   and the timing of Durbin and Koopman.

lh=zeros(1,2);
omega=G*sig*G'+M*M';
% --- robust factorization of the predicted state variance --------------
% omega is symmetric PSD in exact arithmetic; with many exactly observed
% states (MF-BVAR) roundoff makes it slightly asymmetric/indefinite, and
% svd then returns uo ~= vo (or NaN). Symmetrize, use eig, clip negatives.
omega=(omega+omega')/2;
[uo,doo]=eig(omega);
doo=diag(max(real(diag(doo)),0));
vo=uo;
% -------------------------------------------------------------------------
% original codes
%[uo doo vo]=svd(omega);
%--------------------------------------------------------------------------
[u d v]=svd(H*uo*sqrt(doo));
first0=min(find(diag(d)<1e-12));
if isempty(first0),first0=min(size(H))+1;end
u=u(:,1:first0-1);
v=v(:,1:first0-1);
d=diag(d);d=diag(d(1:first0-1));
spred = G*shat; 
fac=vo*sqrt(doo);
yforc = H*spred; 
yhat=y-yforc; 
fhalf=(v/d)*u'; 
fin=fhalf'*fhalf; 
ferr=fhalf*yhat;
lh(1)=-.5*ferr'*ferr;
lh(2)=-sum(log(diag(d)));
kgpart=fac*fhalf; 
% Check 
%comparemat(kgpart,omega*H'*(inv(H*omega*H')) ); 
shatnew=fac*ferr+spred;
signew=fac*(eye(size(v,1))-v*v')*fac';
lh=sum(lh); 