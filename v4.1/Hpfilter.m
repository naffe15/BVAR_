function [s]=hpfilter(y,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function [s]=hpfilter(y,w)
%  Hondrick Prescott filter where:
%  w - smoothing parameter; w=1600 for quarterly data
%  y - the original series that has to be smoothed
%  s - the filtered series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
