% principal components with normalization F'F/T=I
% X is observed
% r is the true number of true factors
% F is T by r matrix of true factors
% Lambda N by r is the true loading matrix
% C=F*Lambda' T by N is the true common component
% chat is the estimated common component

function [ehat,fhat,lambda,ss,Scale]=pc_T(yy,nfac,DEMEAN)

Scale = ones(size(yy,2),1);

if DEMEAN == 2
    [y,Scale]=standard(yy);
elseif DEMEAN ==1 
    y=demean(yy);
else 
    y=yy;
end


[bigt,bign]=size(y);
yy=y*y';
[Fhat0,eigval,Fhat1]=svd(yy);
fhat=Fhat0(:,1:nfac)*sqrt(bigt);
lambda=y'*fhat/bigt;

%chi2=fhat*lambda';
%diag(lambda'*lambda)
%diag(fhat'*fhat)                % this should equal the largest eigenvalues
%sum(diag(eigval(1:nfac,1:nfac)))/sum(diag(eigval))
%mean(var(chi2))                 % this should equal decomposition of variance

ehat=y-fhat*lambda';

ve2=sum(ehat'.*ehat')'/bign;
ss=diag(eigval);

