function bands = compute_rob_cred_bands_(lower_draws, upper_draws, ny, opt)
% COMPUTE_ROB_CRED_BANDS_  Giacomini-Kitagawa robust credible bands.
%
%   bands = compute_rob_cred_bands_(lower_draws, upper_draws, ny, opt)
%
%   Loops over variables (ss = 1:ny), permutes the lower/upper draw arrays
%   into the shape expected by credibleRegion(), and stacks the resulting
%   credible-band matrices into output struct fields .u and .l.
%
%   Inputs
%     lower_draws  -- (ny x hor x ny x K) lower-bound IRF draws
%     upper_draws  -- (ny x hor x ny x K) upper-bound IRF draws
%     ny           -- number of variables
%     opt          -- options struct passed to credibleRegion()
%
%   Output
%     bands.u(:,:,ss)  -- upper credible band for variable ss
%     bands.l(:,:,ss)  -- lower credible band for variable ss

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
