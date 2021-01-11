function [x]=demean(y)
T=size(y,1);
my=repmat(nanmean(y),T,1);
x=(y-my);
