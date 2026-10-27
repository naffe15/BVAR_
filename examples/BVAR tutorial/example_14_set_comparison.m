%% Set identification: old options vs options.set (iresponse_set)
% For each identification scheme the same restrictions are imposed twice:
% with the options of Sections 3.2.3-3.2.6 and with options.set. The script
% reports the computing time and the distance between the IRF bands of the
% identified shocks.
%
% Distance = |percentile(new) - percentile(old)| / width of the old 68%
% band, averaged over variables, horizons and shocks. Draws of
% options.set are summarized with weighted percentiles (wquantile.m) using
% BVAR.irfset_weights. As a benchmark, the script also reports the
% distance due to Monte Carlo noise alone (two halves of the old draws):
% the distance old-new should be of the same order.
%
% Cases: 1) signs; 2) signs + narrative; 3) zeros + signs; 4) long run;
%        5) signs with Giacomini-Kitagawa robust credible sets.
close all; clc; clear;
addpath ../../cmintools/
addpath ../../bvartools/

K           = 1000;   % posterior draws
K_robust    = 200;    % posterior draws for the robust credible sets
L_robust    = 100;    % rotations per draw for the robust credible sets
pp          = [0.16 0.5 0.84];

load DataGK
lags     = 12;
y        = [logip logcpi gs1 ebp];
y_lr     = [diff(logip) logcpi(2:end) gs1(2:end) ebp(2:end)];
varnames = {'IP','CPI','1 year rate','EBP'};

base.hor     = 24;
base.K       = K;
base.noprint = 1;

% ------------------------------------------------------------------ cases
signs = {'y(3,1:3,1)>0','y(2,1:3,1)<0'};          % MP tightening
narr  = {'v([3:5],1)>0','v([11],1)>0'};           % Volcker episodes

C = struct('name',{},'y',{},'old',{},'new',{},'field',{},'shocks',{},'robust',{});

% 1) signs
o = base; o.signs = signs;
n = base; n.set   = signs;
C(end+1) = struct('name','Signs','y',y,'old',o,'new',n, ...
    'field','irsign_draws','shocks',1,'robust',0);

% 2) signs + narrative
o = base; o.signs = signs; o.narrative = narr;
n = base; n.set   = [signs narr];
C(end+1) = struct('name','Signs + narrative','y',y,'old',o,'new',n, ...
    'field','irnarrsign_draws','shocks',1,'robust',0);

% 3) zeros + signs (as in example_3_irf.m, Section 4)
o = base;
o.zeros_signs = {'y(1,1)=1;','y(2,1)=1;','y(3,1)=1;', ...   % AD
                 'y(1,2)=1;','y(2,2)=-1;', ...              % AS
                 'ys(1,3)=0;','ys(2,3)=0;','y(3,3)=1;','y(4,3)=1;'}; % MP
n = base;
n.set = {'y(1,1,1)>0','y(2,1,1)>0','y(3,1,1)>0', ...
         'y(1,1,2)>0','y(2,1,2)<0', ...
         'y(1,1,3)=0','y(2,1,3)=0','y(3,1,3)>0','y(4,1,3)>0'};
C(end+1) = struct('name','Zeros + signs','y',y,'old',o,'new',n, ...
    'field','irzerosign_draws','shocks',1:3,'robust',0);

% 4) long run (Blanchard-Quah: lower triangular long run impact)
o = base; o.long_run_irf = 1;
n = base; n.set = {};
for c = 1 : 4
    for a = 1 : c-1
        n.set{end+1} = sprintf('y(%d,inf,%d)=0', a, c);
    end
    n.set{end+1} = sprintf('y(%d,inf,%d)>0', c, c);
end
C(end+1) = struct('name','Long run','y',y_lr,'old',o,'new',n, ...
    'field','irlr_draws','shocks',1,'robust',0);

% 5) signs with robust credible sets
o = base; o.K = K_robust; o.signs = signs; o.robust_credible_regions.L = L_robust;
n = base; n.K = K_robust; n.set   = signs; n.robust_credible_regions.L = L_robust;
C(end+1) = struct('name','Signs, robust sets','y',y,'old',o,'new',n, ...
    'field','irsign_draws','shocks',1,'robust',1);

% ------------------------------------------------------------------- runs
out = struct('name',{},'t_old',{},'t_new',{},'d_med',{},'d_band',{},'d_noise',{},'ess',{});
for i = 1 : numel(C)
    fprintf('\n%s\n', C(i).name);
    s = C(i).shocks;

    tic; bo = bvar_(C(i).y, lags, C(i).old); t_old = toc;
    tic; bn = bvar_(C(i).y, lags, C(i).new); t_new = toc;
    if strcmp(C(i).name, 'Long run'), b_lr = bo; i_lr = i; end

    Xo = bo.(C(i).field)(:,:,s,:);
    Xn = bn.irfset_draws(:,:,s,:);
    w  = bn.irfset_weights;
    qo = wquantile(Xo, ones(size(Xo,4),1), pp);
    qn = wquantile(Xn, w, pp);
    [d_med, d_band] = band_distance(qo, qn);

    % Monte Carlo noise: distance between the two halves of the old draws,
    % divided by sqrt(2) to match two independent runs with K draws
    % (bvar_ always starts from the same seed, so rerunning it would give
    % identical draws)
    Ko  = size(Xo,4);
    h1  = wquantile(Xo(:,:,:,1:floor(Ko/2)),     ones(floor(Ko/2),1), pp);
    h2  = wquantile(Xo(:,:,:,floor(Ko/2)+1:end), ones(Ko-floor(Ko/2),1), pp);
    [dn1, dn2] = band_distance(h1, h2);
    d_noise = [dn1 dn2] / sqrt(2);
    ess = sum(w)^2 / sum(w.^2);

    out(i) = struct('name',C(i).name,'t_old',t_old,'t_new',t_new, ...
        'd_med',d_med,'d_band',d_band,'d_noise',d_noise,'ess',ess);
    fprintf('  time old %7.1fs  new %7.1fs (x%.1f)\n', t_old, t_new, t_old/t_new);
    fprintf('  distance old-new: median %.3f, 16/84 bands %.3f\n', d_med, d_band);
    fprintf('  MC noise        : median %.3f, 16/84 bands %.3f\n', d_noise(1), d_noise(2));
    fprintf('  effective sample size of the new draws %.0f out of %d\n', ess, sum(w > 0));

    if strcmp(C(i).name, 'Signs + narrative')
        % without the narrative (ADRR) weights, as the old routine
        qu = wquantile(Xn, double(w > 0), pp);
        [dm, db] = band_distance(qo, qu);
        fprintf('  unweighted new draws: median %.3f, 16/84 bands %.3f\n', dm, db);
    end
    if C(i).robust
        ro = bo.irsign_robust_credible_bands;
        rn = bn.irfset_robust_credible_bands;
        wd = ro.u(:,:,s) - ro.l(:,:,s);
        dd = [abs(rn.u(:,:,s)-ro.u(:,:,s)) abs(rn.l(:,:,s)-ro.l(:,:,s))] ./ [wd wd];
        dd = dd(isfinite(dd));
        fprintf('  robust credible sets: distance mean %.3f, max %.3f\n', mean(dd), max(dd));
    end

    % plot: old (grey band, black median) vs new (red dashed)
    figure('Name', C(i).name);
    nv = size(qo,1);
    for sh = 1 : numel(s)
        for v = 1 : nv
            subplot(numel(s), nv, (sh-1)*nv + v); hold on
            hh = 1 : size(qo,2);
            fill([hh fliplr(hh)], [squeeze(qo(v,:,sh,1)) fliplr(squeeze(qo(v,:,sh,3)))], ...
                [.85 .85 .85], 'EdgeColor', 'none');
            plot(hh, squeeze(qo(v,:,sh,2)), 'k', 'LineWidth', 1.5);
            plot(hh, squeeze(qn(v,:,sh,2)), 'r--', 'LineWidth', 1.5);
            plot(hh, squeeze(qn(v,:,sh,1)), 'r:', hh, squeeze(qn(v,:,sh,3)), 'r:');
            title(sprintf('%s, shock %d', varnames{v}, s(sh)));
            axis tight
        end
    end
    legend({'old 68%','old median','new median','new 68%'}, 'Location', 'best');
end

% ---------------------------------------------------------------- summary
fprintf('\n%-22s %9s %9s %8s %10s %10s %10s %10s\n', 'Case', 'old (s)', 'new (s)', 'speedup', ...
    'd median', 'd bands', 'noise med', 'noise bnd');
for i = 1 : numel(out)
    fprintf('%-22s %9.1f %9.1f %8.1f %10.3f %10.3f %10.3f %10.3f\n', out(i).name, out(i).t_old, ...
        out(i).t_new, out(i).t_old/out(i).t_new, out(i).d_med, out(i).d_band, out(i).d_noise);
end

% ------------------------------------------------------- exact check (LR)
% Long run restrictions are point identifying: on the same (Phi,Sigma)
% the two routines must give the same IRF.
Phi   = mean(b_lr.Phi_draws(1:4*lags,:,:),3);
Sigma = mean(b_lr.Sigma_draws,3);
ir_old = iresponse_longrun(Phi, Sigma, base.hor, lags);
ir_new = iresponse_set([], Phi, Sigma, base.hor, C(i_lr).new.set);
fprintf('\nLong run at the posterior mean: max |old - new| = %.2e\n', max(abs(ir_old(:) - ir_new(:))));

% =========================================================================
function [d_med, d_band] = band_distance(qa, qb)
% distance between percentiles, relative to the width of the 68% band of qa
wd = qa(:,:,:,3) - qa(:,:,:,1);
wd(wd < 1e-10) = NaN;                     % e.g. zero restrictions
dm = abs(qb(:,:,:,2) - qa(:,:,:,2)) ./ wd;
db = [abs(qb(:,:,:,1) - qa(:,:,:,1)) ./ wd, abs(qb(:,:,:,3) - qa(:,:,:,3)) ./ wd];
d_med  = mean(dm(isfinite(dm)));
d_band = mean(db(isfinite(db)));
end
