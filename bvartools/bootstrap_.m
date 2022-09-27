function [y,varargout] = bootstrap_(x,Nboot,replacement)

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