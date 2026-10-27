function q = wquantile(X, w, p)
%WQUANTILE Weighted quantiles across draws.
%
%   q = WQUANTILE(X,w,p)
%
%   Inputs
%     X    array of draws, draws on the LAST dimension, e.g.
%          BVAR.irfset_draws (variable x horizon x shock x draws)
%     w    draws x 1   weights, e.g. BVAR.irfset_weights; draws with
%                      w = 0, w = NaN or X = NaN are ignored
%     p    vector      probabilities, e.g. [0.16 0.5 0.84]
%
%   Output
%     q    size of X with the last dimension replaced by numel(p), e.g.
%          variable x horizon x shock x numel(p)
%
%   See also IRESPONSE_SET, PLOT_IRFS_.

%   Filippo Ferroni.

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
