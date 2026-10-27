function [Q,index,flag] = findQs(k,f)
%FINDQS Linear restriction matrices for each shock, from short and long-run zeros.
%
%   [Q,index,flag] = FINDQS(k,f)
%
%   Inputs
%     k       scalar     number of variables
%     f       2k x k     zero-restriction matrix from SIGN2MATRIX (0 marks a
%                        zero restriction)
%
%   Outputs
%     Q       k x 1 cell  restriction matrix for each shock, ordered from
%                         the most to the least restricted
%     index   k x 1       original order of the shocks
%     flag    1 over-identified, 0 exactly identified, -1 under-identified
%
%   See also FINDP, GENERATEDRAW, IRESPONSE_ZEROS_SIGNS.
%
%   Reference: Rubio-Ramirez, J. F., D. F. Waggoner and T. Zha (2010),
%   Review of Economic Studies 77(2), 665-696.

%   Copyright Andrew Binning 2013. Please feel free to use and modify this
%   code as you see fit. If you use it in academic work, please cite
%   Binning, A. (2013), "Underidentified SVAR models: A framework for
%   combining short and long-run restrictions with sign-restrictions",
%   Norges Bank Working Paper 2013/14.

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
