function out = nowcast_bvar(DataMFNC, BVAR, options)
% NOWCAST_BVAR  nowcast_bvar .
%
%   out = nowcast_bvar2(DataMFNC, BVAR, options)
%

fast_kf             = 0;
state_unc           = 0;
KFoptions           = BVAR.KFoptions;
NumbDrawsFromPost   = 1000;
lags                = BVAR.lags;
fhor                = BVAR.fhor;
forecast_data       = BVAR.forecast_data;

[T, ny, nD]         = size(DataMFNC);

if nargin > 2
    if isfield(options,'NumbDrawsFromPost'), NumbDrawsFromPost = options.NumbDrawsFromPost; end
    if isfield(options,'fast_kf'),           fast_kf   = options.fast_kf;                 end
    if isfield(options,'noprint'),           KFoptions.noprint = options.noprint;         end
    if isfield(options,'state_uncertainty'), state_unc = options.state_uncertainty;       end
else
    options = struct();
end

gdp     = KFoptions.mf_varindex;
isflow  = isfield(KFoptions,'index') && numel(KFoptions.index) >= gdp && KFoptions.index(gdp) == 2;
if isflow && lags < 3
    error('nowcast_bvar: the quarterly aggregation of the forecasts needs lags >= 3.');
end

index = randi(BVAR.ndraws, NumbDrawsFromPost, 1);       % same draws as nowcast_bvar for a given seed

NowCast            = nan(T, NumbDrawsFromPost, nD);
frcM_no_shock      = nan(fhor, ny, NumbDrawsFromPost, nD);
frcM_with_shocks   = frcM_no_shock;
frcQ_no_shock      = nan(fhor, NumbDrawsFromPost, nD);
frcQ_with_shocks   = frcQ_no_shock;
Companion_matrix   = diag(ones(ny*(lags-1),1),-ny);

showbar = ~(isfield(KFoptions,'noprint') && KFoptions.noprint == 2);
if showbar, wb = waitbar(0, 'Computing the nowcast with the MF BVAR'); end
KK = NumbDrawsFromPost*nD; dd = 0;

for j = 1:nD
    data = DataMFNC(:,:,j);
    if fast_kf
        tmp = find(sum(isnan(data),2) == 0);
        To  = tmp(end);
    end
    for i = 1:NumbDrawsFromPost
        dd    = dd + 1;
        Phi   = BVAR.Phi_draws(:,:,index(i));
        Sigma = BVAR.Sigma_draws(:,:,index(i));
        Companion_matrix(1:ny,:) = Phi(1:ny*lags,:)';
        opt_i = KFoptions;

        if fast_kf == 1
            % ---- segment 1: rows 1:To (as in nowcast_bvar: BVAR.KFoptions) -----------
            [~, KFzero] = kfilternan(Phi, Sigma, data(1:To,:), KFoptions);
            % ---- segment 2: rows To+1:T, initialised at the end of segment 1 ---------
            opt_i.initialCond = 2;
            opt_i.pZero = KFzero.CovSt(:,:,end);
            opt_i.aZero = KFzero.smoothSt(end,:)';
            [~, KFout] = kfilternan(Phi, Sigma, data(To+1:end,:), opt_i);
            iv  = KFout.index_var(gdp);
            S2  = KFout.smoothSt_plus_ss;                       % (T-To) x ns, levels
            if state_unc
                S2 = S2 + sim_correction(Phi, Sigma, data(To+1:end,:), opt_i, KFout);
            end
            NowCast(1:To,i,j)     = KFzero.smoothSt_plus_ss(:,iv);
            NowCast(To+1:end,i,j) = S2(:,iv);
            Slast = [KFzero.smoothSt_plus_ss; S2];              % full state path, levels
        else
            test = abs(eig(Companion_matrix));
            if any(test > 1.0000000000001), opt_i.initialCond = 1; end   % per draw (not sticky)
            [~, KFout] = kfilternan(Phi, Sigma, data, opt_i);
            iv = KFout.index_var(gdp);
            S2 = KFout.smoothSt_plus_ss;
            if state_unc
                S2 = S2 + sim_correction(Phi, Sigma, data, opt_i, KFout);
            end
            NowCast(:,i,j) = S2(:,iv);
            Slast = S2;
        end

        % ---- forecasts launched from the states at T ------------------------------------
        % first ny columns of the state = current monthly values (latent monthly GDP)
        fd = forecast_data;
        fd.initval = Slast(end-lags+1:end, 1:ny);
        [fns, fws] = forecasts(fd, Phi, Sigma, fhor, lags);
        frcM_no_shock(:,:,i,j)    = fns;
        frcM_with_shocks(:,:,i,j) = fws;
        if isflow
            m0 = fd.initval(end-1:end, gdp);                    % latent monthly GDP at T-1, T
            frcQ_no_shock(:,i,j)    = to_flow(m0, fns(:,gdp));
            frcQ_with_shocks(:,i,j) = to_flow(m0, fws(:,gdp));
        else
            frcQ_no_shock(:,i,j)    = fns(:,gdp);
            frcQ_with_shocks(:,i,j) = fws(:,gdp);
        end
        if showbar, waitbar(dd/KK, wb); end
    end
end
if showbar, close(wb); end

out.NowCast                     = NowCast;
out.options                     = options;
out.state_uncertainty           = state_unc;
out.ForeCast.no_shock           = [NowCast; frcQ_no_shock];
out.ForeCast.with_shocks        = [NowCast; frcQ_with_shocks];
out.ForeCastMonthly.no_shock    = frcM_no_shock;
out.ForeCastMonthly.with_shocks = frcM_with_shocks;
end


% =====================================================================================
function q = to_flow(m0, f)
% quarterly-flow measure (3-month average) of a monthly path:
% m0 = [x(T-1); x(T)], f = [x(T+1); ...; x(T+fhor)]
x = [m0(:); f(:)];
q = (x(3:end) + x(2:end-1) + x(1:end-2)) / 3;
end


% =====================================================================================
function D = sim_correction(Phi, Sigma, data, opt, KFout)
% Durbin-Koopman (2002) simulation smoother, mean-correction form.
% Returns D (T x ns) such that  smoothed mean + D  is a draw of the states
% given the data. Steps: simulate states x+ and observations y+ from the
% state space (same missing pattern as the data), smooth y+, D = x+ - E[x | y+].
[A, B, C, const, Sc] = var2ss(Phi, Sigma, opt.index);   % Sc = chol(Sigma), upper
[T, ny] = size(data);
ns = size(A,1);
% initial state (deviation from steady state)
switch opt.initialCond
    case 2
        a0 = opt.aZero; P0 = opt.pZero;
    case 1
        a0 = zeros(ns,1); P0 = 10*eye(ns);
    otherwise
        a0 = zeros(ns,1); P0 = lyapunov_symm(A, B*(Sc'*Sc)*B');
end
P0 = (P0 + P0')/2;
[V, L] = eig(P0); L = max(real(diag(L)), 0);
x  = a0 + V*(sqrt(L).*randn(ns,1));
% simulate: x(t) = A x(t-1) + B Sc' e(t),  y(t) = C (x(t) + const)
xs = zeros(T, ns);
for t = 1:T
    x = A*x + B*(Sc'*randn(ny,1));
    xs(t,:) = x';
end
ys = (C*(xs' + const))';
ys(isnan(data)) = NaN;
o = opt; o.noprint = 1;
[~, KFs] = kfilternan(Phi, Sigma, ys, o);
D = xs - KFs.smoothSt;                                   % both deviations from steady state
end
