function [K] = fourthmom(x,flag)

[T,N] = size(x);

% %tmp = nan(N*(N+1)/2,N*(N+1)/2);
% tmp_ = nan(N*(N+1)/2,N*(N+1)/2,T);
% for jj = 1 : length(x)
%     tmp = vech(x(jj,:)'*x(jj,:));
%     tmp_(:,:,jj) = tmp * tmp';
%     %tmp_(:,:,jj) = kron(u_(jj,:)'*u_(jj,:),u_(jj,:)'*u_(jj,:));
% end
% K    = mean(tmp_,3);
% 
% end

% default vech; else vec 
if nargin < 2
    flag = 0;
end

%tmp = nan(N*(N+1)/2,N*(N+1)/2);
if flag == 0
tmp_ = nan(N*(N+1)/2,N*(N+1)/2,T);
for jj = 1 : length(x)
    tmp = vech(x(jj,:)'*x(jj,:));
    tmp_(:,:,jj) = tmp * tmp';
    %tmp_(:,:,jj) = kron(u_(jj,:)'*u_(jj,:),u_(jj,:)'*u_(jj,:));
end
elseif flag == 1
    tmp_ = nan(N*N, N*N, T);
    for jj = 1 : length(x)
        tmp = reshape(x(jj,:)'*x(jj,:),N^2,1);
        tmp_(:,:,jj) = tmp * tmp';
    end  
end
K    = mean(tmp_,3);

end
