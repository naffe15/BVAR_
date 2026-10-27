function bands = compute_rob_cred_bands_(lower_draws, upper_draws, ny, opt)
%COMPUTE_ROB_CRED_BANDS_ Robust credible bands for all shocks.
%
%   bands = COMPUTE_ROB_CRED_BANDS_(lower_draws,upper_draws,ny,opt)
%
%   Calls CREDIBLEREGION shock by shock on the lower and upper bounds of
%   the identified sets, e.g. BVAR.irsign_lower_draws and
%   BVAR.irsign_upper_draws.
%
%   Inputs
%     lower_draws  ny x hor x ny x K   lower bounds (variable x horizon x
%                                      shock x draw)
%     upper_draws  ny x hor x ny x K   upper bounds, same layout
%     ny           scalar              number of shocks
%     opt          struct passed to CREDIBLEREGION (.aalpha, .gridLength)
%
%   Output
%     bands.l      ny x hor x ny   lower band (variable x horizon x shock)
%     bands.u      ny x hor x ny   upper band
%
%   See also CREDIBLEREGION, BVAR_.
%
%   Reference: Giacomini, R. and T. Kitagawa (2021), Econometrica 89(4),
%   1519-1556.

bands = struct('u', [], 'l', []);
for ss = 1 : ny
    % flip the order
    % from: variable, horizon, [shock], draw
    % to:   draw    , horizon, [shock], variable
    rMinPost = permute(squeeze(lower_draws(:,:,ss,:)), [3 2 1]);
    rMaxPost = permute(squeeze(upper_draws(:,:,ss,:)), [3 2 1]);
    % Compute robustified credible region (Giacomini-Kitagawa)
    [credlb, credub] = credibleRegion(rMinPost, rMaxPost, opt);
    bands.u(:,:,ss) = credub';
    bands.l(:,:,ss) = credlb';
end
end
