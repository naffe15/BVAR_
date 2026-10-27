function [s]=Hpfilter(y,w)
%HPFILTER Two-sided Hodrick-Prescott filter: returns the trend.
%
%   s = HPFILTER(y,w)
%
%   Inputs
%     y   T x n    series, one per column (a row vector is transposed)
%     w   scalar   smoothing parameter, e.g. 1600 for quarterly data
%
%   Output
%     s   T x n    trend; the cycle is y - s
%
%   See also ONE_SIDED_HPFILTER_SERIAL, HAMFILTER.
%
%   Reference: Hodrick, R. J. and E. C. Prescott (1997), "Postwar U.S.
%   Business Cycles: An Empirical Investigation", Journal of Money, Credit
%   and Banking 29(1), 1-16.

% old function line: function [s]=hpfilter(y,w); renamed to match the file
% name and to avoid a clash with hpfilter of the Econometrics Toolbox

if size(y,1)<size(y,2)
   y=y';
end
t=size(y,1);
a=6*w+1;
b=-4*w;
c=w;
d=[c,b,a];
d=ones(t,1)*d;
m=diag(d(:,3))+diag(d(1:t-1,2),1)+diag(d(1:t-1,2),-1);
m=m+diag(d(1:t-2,1),2)+diag(d(1:t-2,1),-2);
%
m(1,1)=1+w;       m(1,2)=-2*w;
m(2,1)=-2*w;      m(2,2)=5*w+1;
m(t-1,t-1)=5*w+1; m(t-1,t)=-2*w;
m(t,t-1)=-2*w;    m(t,t)=1+w;
%
s=m\y;
