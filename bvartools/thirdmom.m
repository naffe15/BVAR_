function [S] = thirdmom(x,flag)

% default vech; else vec 
if nargin < 2
    flag = 0;
end

[T,N] = size(x);
if flag == 0
    tmp_ = nan(N, N*(N+1)/2, T);
    for jj = 1 : length(x)
        tmp = vech(x(jj,:)'*x(jj,:));
        tmp_(:,:,jj) = x(jj,:)' * tmp';
    end
elseif flag == 1
    tmp_ = nan(N, N*N, T);
    for jj = 1 : length(x)
        tmp = reshape(x(jj,:)'*x(jj,:),N^2,1);
        tmp_(:,:,jj) = x(jj,:)' * tmp';
    end  
end
S    = mean(tmp_,3);

