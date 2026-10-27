function [ehat,fhat,lambda,ss,Scale]=pc_T(yy,nfac,DEMEAN)
%PC_T Principal-components estimate of static factors, normalised F'F/T = I.
%
%   [ehat,fhat,lambda,ss,Scale] = PC_T(yy,nfac)
%   [ehat,fhat,lambda,ss,Scale] = PC_T(yy,nfac,DEMEAN)
%
%   Inputs
%     yy      T x N      panel, time in rows
%     nfac    scalar     number of factors
%     DEMEAN  0, 1 or 2  data treatment: 0 raw (default), 1 demeaned,
%                        2 standardized. Only option 2 fills missing values
%                        (see STANDARD); with 0 or 1 yy must contain no NaN.
%
%   Outputs
%     ehat    T x N      idiosyncratic component, y - fhat*lambda', where y is
%                        the transformed data
%     fhat    T x nfac   factors, fhat'*fhat/T = I
%     lambda  N x nfac   loadings
%     ss      T x 1      eigenvalues of y*y', in decreasing order
%     Scale   N x 1      standard deviations of yy if DEMEAN = 2, ones otherwise
%
%   See also NUMBFACTOR, RESCALEFAVAR, BDFM_.

%   Reference: Bai, J. and S. Ng (2002), "Determining the Number of Factors
%   in Approximate Factor Models", Econometrica 70(1), 191-221.
%   Adapted from the code of Bai and Ng.

Scale = ones(size(yy,2),1);

if nargin < 3 || isempty(DEMEAN)
    DEMEAN = 0;
end
switch DEMEAN
    case 0
        y = yy;
    case 1
        y = demean(yy);
    case 2
        [y,Scale] = standard(yy);
    otherwise
        error('pc_T: DEMEAN must be 0 (raw), 1 (demeaned) or 2 (standardized).');
end
if any(isnan(y(:)))
    error('pc_T: the data contain NaN; use DEMEAN = 2 or fill the missing values first.');
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

