function [x]=demean(y)
%DEMEAN Removes the column means, ignoring NaN.
%
%   x = DEMEAN(y)
%
%   Input
%     y   T x N   data
%
%   Output
%     x   T x N   y minus its column means (nanmean); NaN entries stay NaN
%
%   See also STANDARD.

T=size(y,1);
my=repmat(nanmean(y),T,1);
x=(y-my);
