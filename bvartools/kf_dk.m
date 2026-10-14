function   [shatnew,signew,lh,yhat,fin,kgpart,yforc]=kf_dk(y,H,shat,sig,G,M)

% =========================================================================
% KF_DK  
% function [shatnew,signew,lh,yhat,fin,kgain,yforc]=kf_dk(y,H,shat,sig,G,M)
% 
% This is Chris Sims's KF with a couple of added outputs 
% 1) The (Partial) Kalman Gain and F^(-1) matrices obtained using the 
%     Generalized Inverse are part of the output
% s(t) = G*s(t-1) + R*n(t)  V( n(t) ) = Q 
% then   M=R*Chol(Q)'=R( CQ') 
%        M*M'=R*( CQ'*CQ )*R'
% 
% See KF_MOD. for a related filter 
% NOTE: KGPART is NOT the appropriate Kalman Gain 
% KG=G*KGPART such that 
% KGPART=(P(t)|t-1)*(H')*(F^-1) 
% Use this version when the G matix is time varying and adjust to the 
% timing in DK which have a different timing in the state equation 
% =======================================================================
% Revised, 2/15/2017
% Revised, 3/21/2018

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