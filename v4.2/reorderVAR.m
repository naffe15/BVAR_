function [Phip,Sigmap] = reorderVAR(Phi,Sigma,reordering)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018

% Permute the autoregressive matrix and the variance coveriance of the
% shocks with the new ordering of the variables. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N           = size(Sigma,1);

if length(reordering) ~= N,
    error('Not enough permutations: the vector with new ordering of variables must have size N.')
end

Phip    = Phi(:,reordering);
Sigmap  = nan(N);

for jj  = 1 : length(reordering)
    for hh = 1 : length(reordering)
        
        tmp = Sigma(reordering(jj),reordering(hh));
        Sigmap(jj,hh) = tmp;
    end
    
end

end