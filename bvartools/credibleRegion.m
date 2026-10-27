function [credlb,credub] = credibleRegion(rmin,rmax,opt)
%CREDIBLEREGION Robust credible region for set-identified responses.
%
%   [credlb,credub] = CREDIBLEREGION(rmin,rmax,opt)
%
%   For each horizon and variable, finds the shortest interval that contains
%   the identified set [rmin,rmax] with posterior probability opt.aalpha.
%
%   Inputs
%     rmin    K x H x Ni   posterior draws of the lower bound of the set
%                          (draws x horizon x variable)
%     rmax    K x H x Ni   draws of the upper bound, same layout
%     opt     struct
%       .aalpha       credibility level, e.g. 0.68
%       .gridLength   number of grid points for the centre of the interval
%
%   Outputs
%     credlb  H x Ni   lower end of the robust credible region
%     credub  H x Ni   upper end
%
%   Requires the Statistics Toolbox (quantile).
%
%   See also COMPUTE_ROB_CRED_BANDS_.
%
%   Reference: Giacomini, R. and T. Kitagawa (2021), "Robust Bayesian
%   Inference for Set-Identified Models", Econometrica 89(4), 1519-1556.

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
