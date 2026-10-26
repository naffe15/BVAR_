function q = wquantile(X, w, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'wquantile' computes weighted quantiles across draws, e.g. of the IRFs
% of options.set with the weights BVAR.irfset_weights.
%
% q = wquantile(X, w, p)
%
% Inputs:
% - X, array of draws; the draws are on the LAST dimension
%   (e.g. BVAR.irfset_draws: variable x horizon x shock x draws)
% - w, weights (draws x 1); draws with w = 0 or NaN values are ignored
% - p, probabilities, e.g. [0.16 0.5 0.84]
%
% Output:
% - q, array with the size of X except the last dimension, which has
%   numel(p) elements (e.g. variable x horizon x shock x numel(p))
%
% Filippo Ferroni,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz = size(X);
K  = sz(end);
w  = w(:)';
if numel(w) ~= K
    error('wquantile: the number of weights must equal the number of draws.')
end
Xr = reshape(X, [], K);
q  = nan(size(Xr,1), numel(p));
for i = 1 : size(Xr,1)
    x  = Xr(i,:);
    ok = ~isnan(x) & w > 0 & ~isnan(w);
    if ~any(ok), continue; end
    [xs, ix] = sort(x(ok));
    wo = w(ok);
    cw = cumsum(wo(ix));
    cw = cw / cw(end);
    for j = 1 : numel(p)
        q(i,j) = xs(find(cw >= p(j) - 1e-12, 1));
    end
end
q = reshape(q, [sz(1:end-1) numel(p)]);
end
