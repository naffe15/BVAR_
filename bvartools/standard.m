function [x,Scale]=standard(y)
%STANDARD Standardizes each column to mean 0 and standard deviation 1.
%
%   [x,Scale] = STANDARD(y)
%
%   Input
%     y       T x N   data; NaN values are ignored in the means and std's
%
%   Outputs
%     x       T x N   standardized data; each NaN entry is replaced by an
%                     independent draw from N(0,1)
%     Scale   N x 1   column standard deviations of y
%
%   See also DEMEAN, PC_T.

T=size(y,1);
my=repmat(nanmean(y),T,1);
sy=repmat(nanstd(y),T,1);
x=(y-my)./sy;
inan=isnan(x);
x(inan)=randn(nnz(inan),1);   % one independent N(0,1) draw per missing value
Scale = nanstd(y)';
%x=(y-kron(mean(y),ones(rows(y),1)))./kron(std(y),ones(rows(y),1));
