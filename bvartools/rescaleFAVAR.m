function [ReScale] = rescaleFAVAR(STD,Lambda,n_1,order_pc)

if nargin < 4
    % principal component ordered first
    order_pc = 1;
end

n_w  = size(Lambda,2);
n_2  = size(Lambda,1);

ReScale= NaN;

if size(STD,2)>size(STD,1)
    STD =STD';
end
    
    
switch order_pc
    case 1 % factor first

        Scale_  = repmat([STD; ones(n_1,1)], 1, n_w+n_1);
        Lambda_ = [Lambda zeros(n_2 , n_1); zeros(n_1,n_w+n_1)];
        Lambda_(n_2+1 : n_2+n_1 ,  n_w + 1 : n_w+n_1) = eye(n_1);
        ReScale = Scale_ .* Lambda_ ;

    case 2 % factor second
        
        Scale_  = repmat([ones(n_1,1); STD], 1, n_w +n_1);
        Lambda_ = [zeros(n_1,n_1+n_w); zeros(n_2 , n_1) Lambda];
        Lambda_(1:n_1, 1:n_1) = eye(n_1);
        ReScale = Scale_ .* Lambda_ ;

end