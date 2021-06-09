function [C] = connectedness(Phi,Sigma,p,Omega)
% function [C] = connectedness(Phi,Sigma,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% connectedness 

% Filippo Ferroni, 5/1/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pesaran_Shin = 0;
% dimensions
[m, ny] = size(Phi);    
if rem(m, ny)==0
    lags = m/ny;
else
    lags = floor((m-1)/ny);
end
% if no restriction are imposed use generalized
if nargin<4
    % identification matrix, if nothing is declared use the Pesaran-Shin
    % generalized IRF
    Pesaran_Shin  = 1; 
    Omega = Sigma;        
end

% preaallocation
theta = zeros(ny);
% remove the constant/time trend from the Phi matrix
Phi0 = Phi(1:lags*ny,:);
% compute the MA representation
Psi = var2ma(Phi0,p);
% Equation (1) in Diebold and Yilmaz (2012)
% Construct the Denominator
DD = 0;
for hh = 0 : p-1
    DD = DD + Psi(:,:,hh+1) * Sigma * Psi(:,:,hh+1)';
    % DD = DD + Psi(:,:,hh+1) * Omega * Psi(:,:,hh+1)';
end
% Construct the Numerator and theta
for ss = 1 : ny
    for j= 1 : ny
        NN = 0;
        for hh = 0 : p-1
            %AA = Psi(:,:,hh+1)*Sigma;
            AA = Psi(:,:,hh+1) * Omega;
            NN = NN + AA(ss,j)^2;
        end
        theta(ss,j) = NN / DD(ss,ss);
        if Pesaran_Shin  == 1 
            theta(ss,j) = Sigma(j,j)^(-1/2) * theta(ss,j);            
        end
    end
end

% normalize by the sum of the rows
Theta = theta ./ repmat(sum(theta,2),1,ny);

if max(max(abs(sum(Theta,2) - ones(ny,1))))>1e-8
    error('something went wrong')
end

% remove the main diagonal
% Theta0          = triu(Theta)-tril(Theta);
Theta0          = Theta - diag(diag((Theta)));
C.Index         = sum(sum(Theta0)) / ny * 100; % overall index of connectedness
C.FromAllToUnit = sum(Theta0,2) / (ny-1) * 100;    % connectedness from all to unit i 
C.FromUnitToAll = sum(Theta0,1)' / (ny-1) * 100;   % connectedness from unit i to all 
C.Net           = C.FromUnitToAll - C.FromAllToUnit;   % connectedness from unit i to all 
C.theta         = theta; 
