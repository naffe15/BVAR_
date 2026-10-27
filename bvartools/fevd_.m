function out_ = fevd_(hor,BVAR,options)
%FEVD_ Posterior distribution of the FEVD for every identification in BVAR.
%
%   out_ = FEVD_(hor,BVAR)
%   out_ = FEVD_(hor,BVAR,options)
%
%   Computes FEVD for each posterior draw, with the recursive scheme and,
%   when present in BVAR, with the sign, narrative, zero-sign and
%   heteroskedasticity rotations.
%
%   Inputs
%     hor      scalar   horizon (see FEVD)
%     BVAR     output of BVAR_
%     options  struct (optional)
%       .conf_sig   coverage of the bands (default 0.9)
%
%   Output
%     out_     BVAR with the added fields
%       .fevd             recursive: .mean, .median, .up, .low (ny x ny,
%                         variable x shock, percent) and .horizon
%       .fevdsign, .fevdnarrsign, .fevdzerosign, .fevdheterosked
%                         the same statistics for the other schemes
%
%   See also FEVD, BVAR_.

% alternative identifcation scheme
signs_irf           = 0;
narrative_signs_irf = 0;
zeros_signs_irf     = 0;
heterosked_irf      = 0;
% probability coverage
conf_sig            = 0.9;
% # variables and lags
ny   = BVAR.N;
lags = BVAR.lags;

if nargin >2
    if isfield(options,'conf_sig') == 1
        conf_sig = options.conf_sig;
    end
end
sort_idx   = round((0.5 + [-conf_sig, conf_sig, 0]/2) * BVAR.ndraws);

% recursive FEVD
fevdr = nan(BVAR.N,BVAR.N,BVAR.ndraws);
% actived other FEVD if any
if isempty(BVAR.Omegas) == 0
    signs_irf = 1;
    fevds = nan(BVAR.N,BVAR.N,BVAR.ndraws);
end
if isempty(BVAR.Omegan) == 0
    narrative_signs_irf = 1;
    fevdn = nan(BVAR.N,BVAR.N,BVAR.ndraws);
end
if isempty(BVAR.Omegaz) == 0
    zeros_signs_irf = 1;
    fevdz = nan(BVAR.N,BVAR.N,BVAR.ndraws);
end
if isempty(BVAR.Omegah_draws) == 0
    heterosked_irf = 1;
    fevdh = nan(BVAR.N,BVAR.N,BVAR.ndraws);
end

waitbar_yes = 0;   % old: undefined when BVAR.ndraws <= 99
if BVAR.ndraws > 99
    waitbar_yes = 1;
    wb = waitbar(0, 'Computing the Posterior Distribution of the FEVD');
end
% loop across draws
for nn = 1 : BVAR.ndraws
    Phi   = BVAR.Phi_draws(1 : ny*lags, 1 : ny,nn);
    Sigma = BVAR.Sigma_draws( : , : ,nn);
    % cholesky
    fevdr(:,:,nn) = fevd(hor,Phi,Sigma);
    % sign
    if signs_irf == 1
        fevds(:,:,nn) = fevd(hor,Phi,Sigma,BVAR.Omegas(:,:,nn));
    end
    % sign and narrative
    if narrative_signs_irf == 1
        fevdn(:,:,nn) = fevd(hor,Phi,Sigma,BVAR.Omegan(:,:,nn));
    end
    % zero and sign
    if zeros_signs_irf == 1
        fevdz(:,:,nn) = fevd(hor,Phi,Sigma,BVAR.Omegaz(:,:,nn));
    end
    % heteroskedasticity 
    if heterosked_irf == 1
        % fevdh(:,:,nn) = fevd(hor,Phi,Sigma,BVAR.Omegah(:,:,nn));   % old
        fevdh(:,:,nn) = fevd(hor,Phi,Sigma,BVAR.Omegah_draws(:,:,nn));
    end
    if waitbar_yes, waitbar(nn/BVAR.ndraws, wb); end    
end
if waitbar_yes, close(wb); end

out_ = BVAR;
out_.fevd.horizon = hor; 
% compute statistics and confidence set: recursive
% fevdr_sort = sort(fevdr,4);   % old: sorted along a dimension that does not exist
fevdr_sort = sort(fevdr,3);
out_.fevd.mean   = mean(fevdr,3);
out_.fevd.median = fevdr_sort(:, :, sort_idx(3));
out_.fevd.up     = fevdr_sort(:, :, sort_idx(2));
out_.fevd.low    = fevdr_sort(:, :, sort_idx(1));
% compute statistics and confidence set: sign 
if signs_irf == 1
    out_.fevdsign.mean   = mean(fevds,3);
    % fevds_sort = sort(fevds,4);   % old: sorted along a dimension that does not exist
    fevds_sort = sort(fevds,3);
    out_.fevdsign.median = fevds_sort(:, :, sort_idx(3));
    out_.fevdsign.up     = fevds_sort(:, :, sort_idx(2));
    out_.fevdsign.low    = fevds_sort(:, :, sort_idx(1));
end
% compute statistics and confidence set: sign and narrative
if narrative_signs_irf == 1    
    out_.fevdnarrsign.mean   = mean(fevdn,3);
    % fevdn_sort = sort(fevdn,4);   % old: sorted along a dimension that does not exist
    fevdn_sort = sort(fevdn,3);
    out_.fevdnarrsign.median = fevdn_sort(:, :, sort_idx(3));
    out_.fevdnarrsign.up     = fevdn_sort(:, :, sort_idx(2));
    out_.fevdnarrsign.low    = fevdn_sort(:, :, sort_idx(1));
end
% compute statistics and confidence set: zero and sign
if zeros_signs_irf == 1
    out_.fevdzerosign.mean   = mean(fevdz,3);
    % fevdz_sort = sort(fevdz,4);   % old: sorted along a dimension that does not exist
    fevdz_sort = sort(fevdz,3);
    out_.fevdzerosign.median = fevdz_sort(:, :, sort_idx(3));
    out_.fevdzerosign.up     = fevdz_sort(:, :, sort_idx(2));
    out_.fevdzerosign.low    = fevdz_sort(:, :, sort_idx(1));
end
% compute statistics and confidence set: zero and sign
if heterosked_irf == 1
    out_.fevdheterosked.mean   = mean(fevdh,3);
    % fevdh_sort = sort(fevdh,4);   % old: sorted along a dimension that does not exist
    fevdh_sort = sort(fevdh,3);
    out_.fevdheterosked.median = fevdh_sort(:, :, sort_idx(3));
    out_.fevdheterosked.up     = fevdh_sort(:, :, sort_idx(2));
    out_.fevdheterosked.low    = fevdh_sort(:, :, sort_idx(1));
end


