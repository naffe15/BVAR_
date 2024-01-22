function [credlb,credub] = credibleRegion(rmin,rmax,opt)
% Compute robustified credible region with credibility level alpha, as in
% Kitagawa (2012).
% Inputs:
% - rmin: posterior draws of minimum IR (rows are draws, columns
%         are horizons, third dimension is variable)
% - rmax: posterior draws of maximum IR (same structure as rmin)
% - opt: structure containing options

aalpha = opt.aalpha; % Credibility level
gridLength = opt.gridLength; % Number of points on discrete grid 

[K,H,Ni] = size(rmin);


% Storage.
credlb = zeros(H,Ni);
credub = credlb;

for ii = 1:Ni % For each variable of interest
    
    Cent = zeros(H,1);
    Rad = zeros(H,1);

    for hh = 1:H % For each horizon

        % Construct discrete grid.
        gridr = kron(linspace(min(rmin(:,hh,ii)),max(rmax(:,hh,ii)),...
            gridLength),ones(K,1));

        % d(r,phi) = max{|r-l(phi)|,|r-u(phi)|}.
        d = max(abs(gridr-kron(rmin(:,hh,ii),ones(1,gridLength))),...
            abs(gridr-kron(rmax(:,hh,ii),ones(1,gridLength))));
        zhat = quantile(d,aalpha,1);
        [Rad(hh),ind] = min(zhat,[],2);
        Cent(hh) = gridr(1,ind);

    end
    
    % Approximate robustified credible region is an interval centered at 
    % Cent with radius Rad.
    credlb(:,ii) = Cent - Rad;
    credub(:,ii) = Cent + Rad;

end
    
end
