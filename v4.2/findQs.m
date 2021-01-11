% Copyright Andrew Binning 2013
% Please feel free to use and modify this code as you see if fit. If you
% use this code in any academic work, please cite 
% Andrew Binning, 2013.
% "Underidentified SVAR models: A framework for combining short and long-run restrictions with sign-restrictions,"
% Working Paper 2013/14, Norges Bank.
function [Q,index,flag] = findQs(k,f)
%==========================================================================
% finds the Q matrices that describe the linear restrictions on the shock
% impact matrix. Based on Juan F. Rubio-Ramirez & Daniel F. Waggoner & Tao Zha, 2010.
% "Structural Vector Autoregressions: Theory of Identification and
% Algorithms for Inference," Review of Economic Studies, Oxford University
% Press, vol. 77(2), pages 665-696.
%
% inputs:
% k = number of dependent variables
% f = matrix of short and long run restrictions
%
% outputs:
% Q = a cell that contains the linear restrictions for each equation
% index = the original column order of the matrix of restrictions
% flag = indicates whether the model is over, under or exactly identified
%==========================================================================
E = eye(k);

Q_init = cell(k,2);

for ii = 1:k
    
    Q_init{ii,1} = double(diag(f*E(:,ii)==0));
    Q_init{ii,2} = rank(Q_init{ii,1});
    
end

for ii = 1:k
    
    temp = Q_init{ii,1};
    Q_init{ii,1} = temp(logical(sum(temp,2)),:);
    
end

[new,ord] = sort([Q_init{:,2}],2,'descend');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check identification

if any(new - (k - (1:k)) > 0)  % over-identified
    flag = 1;
elseif all(new - (k - (1:k)) == 0) % exactly identified
    flag = 0;
elseif any(new - (k - (1:k)) < 0) % under-identified
    flag = -1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

index = nan(k,1);

for ii = 1:k
    index(ord(ii)) = ii;
end

Q = cell(k,1);

for ii = 1:k
    
    Q{ii} = Q_init{ord(ii),1};
    
end
