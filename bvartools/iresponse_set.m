function [ir,Omeg,w,info] = iresponse_set(errors,Phi,Sigma,hor,setrestr,cont,opt)
%IRESPONSE_SET Impulse responses of set-identified VARs (options.set).
%
%   [ir,Omeg,w,info] = IRESPONSE_SET(errors,Phi,Sigma,hor,setrestr)
%   [ir,Omeg,w,info] = IRESPONSE_SET(errors,Phi,Sigma,hor,setrestr,cont,opt)
%
%   Any combination of sign, magnitude, zero (short and long run),
%   narrative sign and narrative dominance restrictions, written with the
%   notation of Section 3.2.3 of the Hitchhiker's Guide.
%
%   Inputs
%     errors    T x ny    reduced-form residuals; needed only with narrative
%                         restrictions, otherwise []
%     Phi       (ny*lags+nx) x ny   VAR coefficients (see IRESPONSE)
%     Sigma     ny x ny   covariance of the reduced-form residuals
%     hor       scalar    number of horizons
%     setrestr  cell array of strings (options.set in BVAR_):
%       'y(a,b,c)>m', 'y(a,b,c)<m'  sign (m = 0) or magnitude restriction on
%                                   the response of variable a at horizon b
%                                   to shock c (b = 1 is the impact)
%       'y(a,b,c)=0'                zero restriction at horizon b
%       'y(a,inf,c)=0'              long-run zero restriction
%       'y(a,inf,c)>0'              long-run sign restriction
%                                   a and b can be vectors, e.g.
%                                   'y(3,1:3,1)>0'; c is a scalar
%       any other expression of y   e.g. 'y(1,1,1)/y(2,1,1)>0.5' or
%                                   'max(cumsum(y(1,:,2),2))<1' (finite
%                                   horizons only)
%       'v(tau,n)>0'                narrative sign restriction: shock n is
%                                   positive in periods tau (rows of
%                                   errors); any expression of v, e.g.
%                                   'sum(v(3:5,1))>0'
%       'most(tau,n,a)'             narrative dominance (Type A in
%                                   Antolin-Diaz and Rubio-Ramirez, 2018): shock
%                                   n is the largest contributor to the
%                                   forecast error of variable a in tau
%       'overwhelm(tau,n,a)'        narrative dominance (Type B): shock
%                                   n contributes more than all the others
%                                   together
%     cont      nobs x ny FAVAR only: loadings of the observables; [] for a
%                         standard VAR
%     opt       struct (optional)
%       .draws     accepted rotations to return (default 1)
%       .maxtries  max rotations tried per draw (default 30000)
%       .perm      1 (default): each rotation is checked against every
%                  assignment of its columns and signs to the shocks;
%                  0: sign flips only
%       .M         Monte Carlo draws for the narrative weight (default 5000)
%
%   Outputs
%     ir     ny x hor x ny x draws   responses (variable x horizon x shock x
%                                    draw)
%     Omeg   ny x ny x draws         accepted rotations; the impact matrix is
%                                    chol(Sigma,'lower')*Omeg
%     w      draws x 1               importance weight of each draw (0 when no
%                                    rotation was found); normalise with
%                                    w/sum(w) across posterior draws. It is
%                                    the product of (i) the number of valid
%                                    assignments of columns and signs to the
%                                    shocks, which keeps the prior on the
%                                    rotations uniform, and (ii)
%                                    1/Pr(narrative restrictions hold |
%                                    Phi,Sigma,Omeg) when narrative
%                                    restrictions are given. The Arias,
%                                    Rubio-Ramirez and Waggoner (2018)
%                                    weight for zero restrictions is not
%                                    included.
%     info   struct with fields tries, count (valid assignments), found
%
%   See also WQUANTILE, IRESPONSE_SIGN, IRESPONSE_SIGN_NARRATIVE,
%   IRESPONSE_ZEROS_SIGNS, BVAR_.
%
%   References: Rubio-Ramirez, Waggoner and Zha (2010, REStud); Arias,
%   Rubio-Ramirez and Waggoner (2018, Econometrica); Antolin-Diaz and
%   Rubio-Ramirez (2018, AER), also for the dominance restrictions.

%   Filippo Ferroni.

% ------------------------------------------------------------------ inputs
n = size(Sigma,1);
if nargin < 6 || isempty(cont)
    cont  = eye(n);
    favar = 0;
else
    favar = 1;
end
if nargin < 7 || isempty(opt), opt = struct(); end
if ~isfield(opt,'draws'),    opt.draws    = 1;     end
if ~isfield(opt,'maxtries'), opt.maxtries = 30000; end
if ~isfield(opt,'perm'),     opt.perm     = 1;     end
if ~isfield(opt,'M'),        opt.M        = 5000;  end
if ischar(setrestr), setrestr = {setrestr}; end
if isempty(setrestr)
    error('iresponse_set: no restrictions provided.')
end
nv = size(cont,1);
T  = size(errors,1);
doweight = (nargout > 2);

% ------------------------------------------------------- parse restrictions
% (cached: the same restrictions are parsed once across posterior draws)
persistent cache_key cache_R
key = [sprintf('%s|', setrestr{:}) sprintf('%d,', [n nv hor T])];
if isequal(key, cache_key)
    R = cache_R;
else
    R = parse_set(setrestr, n, nv, hor, T);
    cache_key = key;
    cache_R   = R;
end
if (~isempty(R.N) || ~isempty(R.D)) && T == 0
    error('iresponse_set: narrative restrictions need the reduced form errors.')
end

% ------------------------------------- restrictions as matrices on Omega
% y(:,h,:) = IR(:,:,h) * Omega, with IR(:,:,h) = cont * Psi_{h-1} * A
A     = cholA(Sigma);
hneed = max([1, R.lin.b(isfinite(R.lin.b)), R.hgen]);
ir0   = iresponse(Phi,Sigma,hneed,eye(n));
IR    = zeros(nv,n,hneed);
for h = 1 : hneed
    IR(:,:,h) = cont * reshape(ir0(:,h,:),n,n);
end
if any(isinf(R.lin.b))
    IRinf = cont * long_run(Phi,A);
end
Z  = cell(1,n);     % zero restrictions:       Z{c} * q_c = 0
Sd = cell(1,n);     % inequality restrictions: Sd{c} * q_c > md{c}
md = cell(1,n);
for r = 1 : numel(R.lin.a)
    if isinf(R.lin.b(r)), row = IRinf(R.lin.a(r),:);
    else,                 row = IR(R.lin.a(r),:,R.lin.b(r));
    end
    c = R.lin.c(r);
    if R.lin.zero(r)
        Z{c}  = [Z{c}; row];
    else
        Sd{c} = [Sd{c}; R.lin.d(r)*row];
        md{c} = [md{c}; R.lin.d(r)*R.lin.m(r)];
    end
end

% narrative restrictions 'v(tau,n) OP m' are linear too: v = ustd * Omega
if ~isempty(R.N) || ~isempty(R.D)
    ustd = errors / A';
end
for i = find([R.N.lin])
    c = R.N(i).shock1;
    for t = R.N(i).rows
        Sd{c} = [Sd{c}; R.N(i).d * ustd(t,:)];
        md{c} = [md{c}; R.N(i).d * R.N(i).m];
    end
end

% order in which the columns are drawn (most zeros first) and rank check
nz        = cellfun(@(x) size(x,1), Z);
[nzs,ord] = sort(nz,'descend');
if any(nzs > n - (1:n))
    error('iresponse_set: too many zero restrictions (rank condition fails).')
end
haszeros = any(nz > 0);

% shocks checked column by column (single-shock restrictions)
Nl   = R.N([R.N.lin]);                       % linear: checked through Sd
Ns   = R.N([R.N.single] & ~[R.N.lin]);         % one shock, generic
Nx   = R.N(~[R.N.single]);                     % several shocks
Lset = unique([R.lin.c(~R.lin.zero), [R.gy.shock], [Ns.shock1], [Nl.shock1], [R.D.shock]]);
pool = find(nz == 0);                  % columns that can be permuted
Lp   = Lset(nz(Lset) == 0);            % checked shocks drawn from the pool
Lz   = Lset(nz(Lset) >  0);            % checked shocks with zeros (own column)
k    = numel(Lp);
p    = numel(pool);
if k > p
    error('iresponse_set: more identified shocks than free columns.')
end
% all injections of the k pool shocks into the p pool columns
if opt.perm && k > 0
    np = prod(p-k+1 : p);
    if np > 1e6
        warning('iresponse_set: too many permutations, only sign flips are used.')
        opt.perm = 0;
    end
end
if opt.perm && k > 0
    cmb   = nchoosek(1:p, k);
    Plist = zeros(size(cmb,1)*factorial(k), k);
    row   = 0;
    for i = 1 : size(cmb,1)
        P = perms(cmb(i,:));
        Plist(row+1:row+size(P,1),:) = P;
        row = row + size(P,1);
    end
else
    [~,Plist] = ismember(Lp, pool);    % shock c keeps column c
end
np    = size(Plist,1);
if k > 0
    lin_k = sub2ind([k p], repmat(1:k,np,1), Plist);
else
    Plist = zeros(1,0);
    np    = 1;
end
gy_sh = [R.gy.shock];
Ns_sh = [Ns.shock1];
D_sh  = [R.D.shock];

% ----------------------------------------------------------------- outputs
nd    = opt.draws;
ir    = nan(nv,hor,n,nd);
Omeg  = nan(n,n,nd);
w     = zeros(nd,1);
info  = struct('found',zeros(nd,1),'tries',opt.maxtries*ones(nd,1),'count',zeros(nd,1));

% -------------------------------------------------------------- main loop
for dd = 1 : nd
    for tr = 1 : opt.maxtries
        % ---- candidate rotation (zeros imposed column by column)
        if haszeros
            Q = draw_zeros(Z, ord, n);
        else
            Q = generateQ(n);
        end

        % ---- check every candidate column of every checked shock
        OP   = false(n,n);             % OP(c,j): column  j fits shock c
        ON   = false(n,n);             % ON(c,j): column -j fits shock c
        yQ   = [];  vall = [];  Bo = [];
        fail = false;
        for c = Lset
            if nz(c) > 0, cols = c; else, cols = pool; end
            okP = true(1,numel(cols));
            okN = okP;
            if ~isempty(Sd{c})
                X   = Sd{c} * Q(:,cols);
                okP = all(X  > md{c}, 1);
                okN = all(-X > md{c}, 1);
            end
            % generic expressions of y involving only shock c
            ig = find(gy_sh == c);
            if ~isempty(ig) && any(okP | okN)
                if isempty(yQ), yQ = irf_rot(IR, Q, R.hgen); end
                for jj = find(okP | okN)
                    y = zeros(nv,R.hgen,n);
                    if okP(jj)
                        y(:,:,c) = yQ(:,:,cols(jj));
                        okP(jj)  = eval_y(R.gy(ig), y);
                    end
                    if okN(jj)
                        y(:,:,c) = -yQ(:,:,cols(jj));
                        okN(jj)  = eval_y(R.gy(ig), y);
                    end
                end
            end
            % narrative sign restrictions involving only shock c
            in = find(Ns_sh == c);
            if ~isempty(in) && any(okP | okN)
                if isempty(vall), vall = ustd * Q; end
                V  = zeros(T,n);
                sm = ones(1,n); sm(c) = -1;
                for jj = find(okP | okN)
                    V(:,c) = vall(:,cols(jj));
                    if okP(jj), okP(jj) = check_narrative(Ns(in), R.D([]), V, (1:T)', ones(1,n), []); end
                    if okN(jj), okN(jj) = check_narrative(Ns(in), R.D([]), V, (1:T)', sm, []); end
                end
            end
            % narrative dominance (depends only on the column, not its sign)
            id = find(D_sh == c);
            if ~isempty(id) && any(okP | okN)
                if isempty(vall), vall = ustd * Q; end
                if isempty(Bo),   Bo   = cont * A * Q; end
                okD = true(1,n);
                for i = id
                    for t = R.D(i).tau
                        cc = abs(Bo(R.D(i).var,:) .* vall(t,:));
                        if R.D(i).type == 'A'
                            [m1,i1] = max(cc);
                            c2 = cc; c2(i1) = -inf;
                            okt = false(1,n);
                            okt(i1) = m1 > max(c2);
                        else
                            okt = cc > sum(cc) - cc;
                        end
                        okD = okD & okt;
                    end
                end
                okP = okP & okD(cols);
                okN = okN & okD(cols);
            end
            if ~any(okP | okN), fail = true; break; end
            OP(c,cols) = okP;
            ON(c,cols) = okN;
        end
        if fail, continue; end
        Vn = double(OP) + double(ON);  % number of valid signs

        % ---- valid assignments of pool columns to the pool shocks
        if k > 0
            VV    = Vn(Lp, pool);
            prodV = prod(reshape(VV(lin_k), np, k), 2);
        else
            prodV = 1;
        end
        cntz  = prod(Vn(sub2ind([n n], Lz, Lz)));
        total = sum(prodV) * cntz;
        if total == 0, continue; end
        pick  = find(cumsum(prodV) >= rand*sum(prodV), 1);

        % ---- build the rotation
        Qn = Q;
        for c = Lz
            Qn(:,c) = pick_sign(OP(c,c), ON(c,c)) * Q(:,c);
        end
        used = pool(Plist(pick,:));
        if k == 0, used = []; end
        for i = 1 : k
            Qn(:,Lp(i)) = pick_sign(OP(Lp(i),used(i)), ON(Lp(i),used(i))) * Q(:,used(i));
        end
        free_sh  = setdiff(pool, Lp);           % unrestricted shocks
        free_col = setdiff(pool, used);
        Qn(:,free_sh) = Q(:,free_col);

        % ---- restrictions involving several shocks
        if ~isempty(R.gx)
            y = irf_rot(IR, Qn, R.hgen);
            if ~eval_y(R.gx, y), continue; end
        end
        if ~isempty(Nx)
            V = ustd * Qn;
            if ~check_narrative(Nx, R.D([]), V, (1:T)', ones(1,n), []), continue; end
        end

        % ---- accepted
        Omeg(:,:,dd) = Qn;
        y0 = iresponse(Phi,Sigma,hor,Qn);
        if favar
            for jj = 1 : n, y0(:,:,jj) = cont * y0(:,:,jj); end
            ir(:,:,:,dd) = y0;
        else
            ir(:,:,:,dd) = y0;
        end
        info.found(dd) = 1;
        info.tries(dd) = tr;
        info.count(dd) = total;
        w(dd) = total;
        if doweight && (~isempty(R.N) || ~isempty(R.D))
            w(dd) = total * narrative_weight(R.N, R.D, cont*A*Qn, n, T, opt.M);
        end
        break
    end
end

if any(info.found == 0)
    warning('I could not find a rotation satisfying the restrictions.')
end

end

% =========================================================================
%                              sub-functions
% =========================================================================

function s = pick_sign(okp, okn)
if okp && okn
    s = 1 - 2*(rand < 0.5);
elseif okp
    s = 1;
else
    s = -1;
end
end

% -------------------------------------------------------------------------
function Q = draw_zeros(Z, ord, n)
% Columns drawn one at a time, uniformly on the unit sphere of the space
% orthogonal to the zero restrictions and to the columns already drawn
% (Rubio-Ramirez, Waggoner and Zha 2010; Arias et al. 2018).
Q = zeros(n);
for i = 1 : n
    c  = ord(i);
    Mc = [Z{c}; Q(:,ord(1:i-1))'];
    if isempty(Mc)
        Nc = eye(n);
    else
        Nc = null(Mc);
    end
    x      = randn(size(Nc,2),1);
    Q(:,c) = Nc * (x / norm(x));
end
end

% -------------------------------------------------------------------------
function y = irf_rot(IR, Q, hg)
% y(:,h,:) = IR(:,:,h) * Q for h = 1,...,hg
nv = size(IR,1);
n  = size(Q,2);
y  = zeros(nv,hg,n);
for h = 1 : hg
    y(:,h,:) = reshape(IR(:,:,h) * Q, nv, 1, n);
end
end

% -------------------------------------------------------------------------
function L = long_run(Phi, A)
% long run impact (I - A_1 - ... - A_p)^(-1) * A
n = size(A,1);
m = size(Phi,1);
if rem(m,n) == 0, lags = m/n; else, lags = floor((m-1)/n); end
S = zeros(n);
for l = 1 : lags
    S = S + Phi((l-1)*n+1 : l*n, 1:n)';
end
L = (eye(n) - S) \ A;
end

% -------------------------------------------------------------------------
function A = cholA(Sigma)
% same factorization as iresponse.m
try
    A = chol(Sigma,'lower');
catch
    iS        = pinv(Sigma);
    [L, D, P] = ldl(iS);
    A         = (sqrt(D) * L' * P')';
end
end

% -------------------------------------------------------------------------
function pass = eval_y(S, y) %#ok<INUSD> y is used inside eval
pass = true;
for k = 1 : numel(S)
    tmp = eval(S(k).str);
    if isempty(tmp) || any(isnan(tmp(:)))
        error('There is a nan in the restriction %s.', S(k).str);
    end
    if ~all(tmp(:)), pass = false; return; end
end
end

% -------------------------------------------------------------------------
function ok = check_narrative(N, D, V, map, s, Bo)
% Narrative sign and dominance restrictions.
% V   : (#rows) x n x M structural shocks (M = 1: realized shocks;
%       M > 1: Monte Carlo draws for the importance weight)
% map : T x 1, maps the period t into the row of V
% s   : 1 x n sign of each shock
M  = size(V,3);
n  = size(V,2);
ok = true(1,M);
v  = @(t,j) vget(V, map, s, t, j); %#ok<NASGU> v is used inside eval
for k = 1 : numel(N)
    tmp = eval(N(k).str);
    if isempty(tmp) || any(isnan(tmp(:)))
        error('There is a nan in the narrative restriction %s.', N(k).str);
    end
    if M == 1
        ok = ok & all(tmp(:));
    else
        ok = ok & all(reshape(tmp,[],M),1);
    end
    if ~any(ok), return; end
end
for k = 1 : numel(D)
    a  = D(k).var;
    j  = D(k).shock;
    for t = D(k).tau
        Ev = reshape(V(map(t),:,:), n, M);
        cc = abs(Bo(a,:).' .* Ev);
        cj = cc(j,:);
        cc(j,:) = [];
        if isempty(cc), continue; end
        if D(k).type == 'A'
            ok = ok & (cj > max(cc,[],1));
        else
            ok = ok & (cj > sum(cc,1));
        end
    end
    if ~any(ok), return; end
end
end

% -------------------------------------------------------------------------
function out = vget(V, map, s, t, j)
% v(t,j) inside the narrative restrictions: (numel(t)*numel(j)) x M
if ischar(t), t = 1 : numel(map); end
if ischar(j), j = 1 : size(V,2);  end
r   = map(t);
X   = V(r, j, :) .* reshape(s(j), 1, []);
out = reshape(X, numel(r)*numel(j), size(V,3));
end

% -------------------------------------------------------------------------
function w = narrative_weight(N, D, Bo, n, T, M)
% Antolin-Diaz and Rubio-Ramirez (2018): 1/Pr(narrative restrictions hold),
% probability over v_t ~ N(0,I) on the narrative periods. Restrictions
% that share no period are independent: the probability is computed for
% each block and multiplied.
nN   = numel(N);
rows = [{N.rows}, {D.tau}];
nr   = numel(rows);
lab  = 1 : nr;
changed = true;
while changed
    changed = false;
    for i = 1 : nr
        for k = i+1 : nr
            if lab(i) ~= lab(k) && ~isempty(intersect(rows{i}, rows{k}))
                lab(lab == max(lab(i),lab(k))) = min(lab(i),lab(k));
                changed = true;
            end
        end
    end
end
blocks = unique(lab);
prob   = zeros(1,numel(blocks));
for b = 1 : numel(blocks)
    in  = find(lab == blocks(b));
    iN  = in(in <= nN);
    iD  = in(in >  nN) - nN;
    rr  = unique([rows{in}]);
    map = zeros(T,1);
    map(rr) = 1 : numel(rr);
    V   = randn(numel(rr), n, M);
    okm = check_narrative(N(iN), D(iD), V, map, ones(1,n), Bo);
    prob(b) = max(mean(okm), 1/M);
end
w = 1 / prod(prob);
end

% =========================================================================
%                                 parser
% =========================================================================

function R = parse_set(str, n, nv, hor, T)
R.lin  = struct('a',[],'b',[],'c',[],'zero',false(1,0),'d',[],'m',[]);
R.gy   = struct('str',{},'shock',{});      % generic, one shock
R.gx   = struct('str',{});                 % generic, several shocks
R.N    = struct('str',{},'rows',{},'shocks',{},'single',{},'shock1',{},'lin',{},'d',{},'m',{});
R.D    = struct('str',{},'tau',{},'shock',{},'var',{},'type',{});
R.hgen = 1;
for k = 1 : numel(str)
    s = strtrim(str{k});
    if ~isempty(regexp(s, '^(most|overwhelm)\s*\(', 'once'))
        R.D(end+1) = parse_dominant(s, n, T, nv);
    elseif ~isempty(find_calls(s,'v'))
        if ~isempty(find_calls(s,'y'))
            error('Restriction ''%s'': y and v cannot be mixed.', s);
        end
        R.N(end+1) = parse_narrative(s, n, T);
    elseif ~isempty(find_calls(s,'y'))
        [isl, a, b, c, op, m] = parse_linear(s, n, nv, hor);
        if isl
            if any(strcmp(op, {'=','=='}))
                if m ~= 0
                    error('Restriction ''%s'': only zero equality restrictions are allowed.', s);
                end
                zr = true;  d = 1;
            else
                zr = false;
                d  = 1 - 2*any(strcmp(op, {'<','<='}));
            end
            [aa, bb] = ndgrid(a, b);
            nr = numel(aa);
            R.lin.a    = [R.lin.a,    aa(:)'];
            R.lin.b    = [R.lin.b,    bb(:)'];
            R.lin.c    = [R.lin.c,    c*ones(1,nr)];
            R.lin.zero = [R.lin.zero, repmat(zr,1,nr)];
            R.lin.d    = [R.lin.d,    d*ones(1,nr)];
            R.lin.m    = [R.lin.m,    m*ones(1,nr)];
        else
            if ~isempty(regexp(s, '(?<![<>=~])=(?!=)', 'once'))
                error(['Restriction ''%s'': zero restrictions must be written ' ...
                    'as ''y(a,b,c)=0''.'], s);
            end
            calls  = find_calls(s,'y');
            shocks = [];
            for i = 1 : numel(calls)
                ar = calls{i};
                if numel(ar) ~= 3
                    error('Restriction ''%s'': use y(a,b,c).', s);
                end
                cc = eval_index(ar{3}, n);
                hh = eval_index(ar{2}, hor);
                if isempty(cc) || isempty(hh) || any(~isfinite(hh))
                    error(['Restriction ''%s'': could not read the indexes ' ...
                        '(inf is only allowed in y(a,inf,c)>m, <m or =0).'], s);
                end
                shocks = [shocks cc(:)']; %#ok<AGROW>
                R.hgen = max(R.hgen, max(hh));
            end
            shocks = unique(shocks);
            if numel(shocks) == 1
                R.gy(end+1) = struct('str', s, 'shock', shocks);
            else
                R.gx(end+1) = struct('str', s);
            end
        end
    else
        error('Restriction ''%s'': it must contain y(...), v(...), most(...) or overwhelm(...).', s);
    end
end
end

% -------------------------------------------------------------------------
function [isl, a, b, c, op, m] = parse_linear(s, n, nv, hor)
% 'y(a,b,c) OP m' with a single call to y on the left and a number on
% the right
isl = false; a = []; b = []; c = []; op = ''; m = [];
if ~strncmp(s, 'y', 1), return; end
[calls, ends] = find_calls(s, 'y');
if numel(calls) ~= 1 || numel(calls{1}) ~= 3, return; end
tok = regexp(s(ends(1)+1:end), '^\s*(>=|<=|==|=|>|<)\s*(.+?)\s*;?\s*$', 'tokens', 'once');
if isempty(tok), return; end
mm = str2double(tok{2});
if isnan(mm)
    try, mm = eval(tok{2}); catch, return; end
end
if ~isnumeric(mm) || numel(mm) ~= 1 || ~isfinite(mm), return; end
ar = calls{1};
a  = eval_index(ar{1}, nv);
b  = eval_index(ar{2}, hor);
c  = eval_index(ar{3}, n);
if isempty(a) || isempty(b) || numel(c) ~= 1
    error('Restriction ''%s'': could not read the indexes.', s);
end
if any(a < 1 | a > nv) || c < 1 || c > n || any(b < 1)
    error('Restriction ''%s'': index out of range.', s);
end
isl = true; op = tok{1}; m = mm;
end

% -------------------------------------------------------------------------
function N = parse_narrative(s, n, T)
s0    = regexprep(s, '\<end\>', num2str(T));
calls = find_calls(s0, 'v');
rows = []; shocks = [];
for i = 1 : numel(calls)
    a = calls{i};
    if numel(a) ~= 2
        error(['Narrative restriction ''%s'': use v(tau,n), e.g. ' ...
            '''sum(v(3:5,1))>0''.'], s);
    end
    t = eval_index(a{1}, T);
    c = eval_index(a{2}, n);
    if isempty(t) || isempty(c) || any(t < 1 | t > T | t ~= round(t))
        error('Narrative restriction ''%s'': periods must be between 1 and %d.', s, T);
    end
    rows   = [rows t(:)'];   %#ok<AGROW>
    shocks = [shocks c(:)']; %#ok<AGROW>
end
shocks = unique(shocks);
N = struct('str', s0, 'rows', unique(rows), 'shocks', shocks, ...
    'single', numel(shocks) == 1, 'shock1', shocks(1), 'lin', false, 'd', 0, 'm', 0);
% linear form 'v(tau,n) OP m' (one call to v, a number on the right)
[ca, en] = find_calls(s0, 'v');
if numel(ca) == 1 && strncmp(strtrim(s0), 'v', 1)
    tok = regexp(s0(en(1)+1:end), '^\s*(>=|<=|>|<)\s*(.+?)\s*;?\s*$', 'tokens', 'once');
    if ~isempty(tok) && numel(N.shocks) == 1
        mm = str2double(tok{2});
        if isnan(mm)
            try, mm = eval(tok{2}); catch, mm = NaN; end
        end
        if isnumeric(mm) && numel(mm) == 1 && isfinite(mm)
            N.lin = true;
            N.d   = 1 - 2*any(strcmp(tok{1}, {'<','<='}));
            N.m   = mm;
        end
    end
end
end

% -------------------------------------------------------------------------
function D = parse_dominant(s, n, T, nv)
s0    = regexprep(s, '\<end\>', num2str(T));
tok   = regexp(s0, '^(most|overwhelm)', 'tokens', 'once');
calls = find_calls(s0, tok{1});
a     = calls{1};
if numel(a) ~= 3
    error('Dominance restriction ''%s'': use most(tau,n,a) or overwhelm(tau,n,a).', s);
end
t = eval_index(a{1}, T);
c = eval_index(a{2}, n);
v = eval_index(a{3}, nv);
if isempty(t) || numel(c) ~= 1 || numel(v) ~= 1 || any(t < 1 | t > T)
    error('Dominance restriction ''%s'': could not read the indexes.', s);
end
if strcmp(tok{1},'most'), ty = 'A'; else, ty = 'B'; end
D = struct('str', s, 'tau', unique(t(:)'), 'shock', c, 'var', v, 'type', ty);
end

% -------------------------------------------------------------------------
function x = eval_index(s, nmax)
s = strtrim(s);
if strcmp(s, ':'), x = 1 : nmax; return; end
try
    x = eval(s);
catch
    x = [];
end
end

% -------------------------------------------------------------------------
function [args, ends] = find_calls(str, name)
% every occurrence of name(...) in str: top-level arguments as strings
% and position of the closing parenthesis
args = {};
ends = [];
idx  = regexp(str, ['(?<![\w\.])' name '\s*\(']);
for i = idx
    p     = i + numel(name);
    while str(p) ~= '(', p = p + 1; end
    depth = 0;
    start = p + 1;
    parts = {};
    for q = p : numel(str)
        ch = str(q);
        if any(ch == '([{')
            depth = depth + 1;
        elseif any(ch == ')]}')
            depth = depth - 1;
            if depth == 0
                parts{end+1} = str(start:q-1); %#ok<AGROW>
                break
            end
        elseif ch == ',' && depth == 1
            parts{end+1} = str(start:q-1); %#ok<AGROW>
            start = q + 1;
        end
    end
    args{end+1} = parts; %#ok<AGROW>
    ends(end+1) = q;     %#ok<AGROW>
end
end
