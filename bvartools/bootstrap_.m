function [y,varargout] = bootstrap_(x,Nboot,replacement)
%BOOTSTRAP_ Resamples the rows of a matrix.
%
%   y = BOOTSTRAP_(x,Nboot)
%   [y,order_] = BOOTSTRAP_(x,Nboot,replacement)
%
%   Inputs
%     x            T x N   data, e.g. residuals
%     Nboot        scalar  number of samples
%     replacement  0 (default) random permutation of the rows, 1 draw
%                  with replacement
%
%   Outputs
%     y        T x N x Nboot   resampled data
%     order_   T x Nboot       rows used in each sample
%
%   See also CVAR_.

if nargin<3
    replacement = 0;
else
    if replacement>1
        error('replacement is logical, 0 or 1.')
    end
end

[T,N] = size(x);

y = nan(T,N,Nboot);

for nboot =1 :Nboot
    if replacement ==1
        % replacement
        order_(:,nboot)  = randi(T,[T,1]);        
    else 
        % no replacement
        order_(:,nboot)  = randperm(T);
    end
    y(:,:,nboot) = x(order_(:,nboot),:);
end


varargout = {order_};