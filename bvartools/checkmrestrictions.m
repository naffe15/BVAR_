function d = checkmrestrictions(hmoments,v,f)

d = 0;
% set of higher moments - contains the moments of 
% 1st dimension:   shocks
% 2nd dimension:   moment estimator (1= based on sample 3rd and 4th
% moments, #>1 robust estimators
% 3rd dimension:   moment order (1 = skewnewss, 2 = kurtosis)
nv = size(v,2);
hm = nan(nv,4,3);

% excess skewness - 
[skw,skw_normal]  = skewness_(v,f(:,:,1));
hm(:,:,1) = skw - repmat(skw_normal',nv,1);

% excess kurtosis - 
[krt,krt_normal]  = kurtosis_(v,f(:,:,2));        
hm(:,:,2) = krt - repmat(krt_normal,nv,1);

% % right tails only - 
% [krt,krt_normal]  = kurtosis_(v,f(:,:,2));        
% hm(:,:,3) = krt - repmat(krt_normal,nv,1);

% excess LEFT tail kurtosis - 
% TBA

% excess RIGHT tail kurtosis - 
% TBA

d=0; 
count   = 0;
for ii = 1 : size(hmoments,2)
    tmp = eval(hmoments{ii});
    count = count + min(tmp);
end
if isempty(count)==1  
    error('There is a nan in the restrictions.');
end
if count == size(hmoments,2) % if all signs are verified stop
    d=1;
end

