function result = gda_structural_scenario_moments(Phi,Sigma,Q,forecast_data,scenario)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'gda_structural_scenario_moments' computes conditional moments and draws for
% one SVAR scenario and one parameter/identification draw
% Called by bvar_ inside its draw loop when options.scenario is set (%%gda blocks).
%
% Reference: Antolin-Diaz, J., I. Petrella and J. F. Rubio-Ramirez (2021), "Structural
% scenario analysis with SVARs", Journal of Monetary Economics 117, 798-815.
%
% Inputs:
% - Phi, reduced-form VAR parameters
% - Sigma, reduced-form covariance matrix
% - Q, orthogonal identification matrix
% - forecast_data, initial values and future exogenous data
% - scenario, structure with horizon, C, f, driving_shocks, uncertainty,
%   Omega_f, n_draws and name
%
% Output:
% - result, conditional moments, paths, restrictions and plausibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Scenario setup
% Read dimensions and defaults used by later blocks.
fhor = scenario.horizon;
Nvar = size(Sigma,1);
Nhor = Nvar*fhor;
lags = size(forecast_data.initval,1);
if norm(Q'*Q-eye(Nvar),'fro') > 1e-8*Nvar
    error('scenario:rotation','Q must be an orthogonal rotation.');
end

if ~isfield(scenario,'C') || isempty(scenario.C)
    C = zeros(0,Nhor);
else
    C = scenario.C;
end
if ~isfield(scenario,'f')
    f = zeros(0,1);
else
    f = scenario.f(:);
end

if ~isfield(scenario,'driving_shocks')
    driving = 1:Nvar;
else
    driving = unique(scenario.driving_shocks(:)','stable');
end
non_driving = setdiff(1:Nvar,driving,'stable');

if ~isfield(scenario,'n_draws')
    Ndraws = 0;
else
    Ndraws = scenario.n_draws;
end

%% 2. No-shock forecast
% Compute b, the deterministic component of the forecast (Antolin-Diaz et al. 2021, equation (3), definition in Appendix A).
lags_data = forecast_data.initval(:,:,1);
sims_no_shock = zeros(fhor,Nvar);
for t = 1:fhor
    X = [reshape(flip(lags_data,1)',1,Nvar*lags) forecast_data.xdata(t,:)];
    y = X*Phi;
    sims_no_shock(t,:) = y;
    lags_data = [lags_data(2:end,:); y];
end
b = reshape(sims_no_shock',Nhor,1);

%% 3. Structural transmission matrix
% Build the impulse responses and G, the matrix that maps the stacked structural shocks into the future path
% (Antolin-Diaz et al. 2021, equation (3), definition in Appendix A).
impact = chol((Sigma+Sigma')/2,'lower')*Q;
B = reshape(Phi(1:Nvar*lags,:)',Nvar,Nvar,lags);
Psi = zeros(Nvar,Nvar,fhor);
Psi(:,:,1) = eye(Nvar);
for t = 2:fhor
    for lag = 1:min(lags,t-1)
        Psi(:,:,t) = Psi(:,:,t)+B(:,:,lag)*Psi(:,:,t-lag);
    end
end

G = zeros(Nhor);
for t = 1:fhor
    rows = (t-1)*Nvar+(1:Nvar);
    for j = 1:t
        cols = (j-1)*Nvar+(1:Nvar);
        G(rows,cols) = Psi(:,:,t-j+1)*impact;
    end
end

%% 4. Uncertainty around observable restrictions
% Select Omega_f for the hard, variance-preserving or user-defined case.
Ncond = size(C,1);
if Ncond == 0
    Omega_f = zeros(0);
else
    policy = 'user';
    if isfield(scenario,'uncertainty') && ~isempty(scenario.uncertainty)
        policy = lower(char(scenario.uncertainty));
    end
    switch policy
        case 'hard'
            Omega_f = zeros(Ncond);
        case 'preserve_unconditional'
            Omega_f = C*(G*G')*C';
        case 'user'
            Omega_f = scenario.Omega_f;
    end
    Omega_f = (Omega_f+Omega_f')/2;
end

%% 5. Structural scenario moments
% Stack the restrictions on the observables and on the non-driving shocks (Antolin-Diaz et al. 2021, equation (15))
% and apply the Penrose (minimum norm) solution (equations (10)-(14)).
S = eye(Nvar);
Xi = kron(eye(fhor),S(non_driving,:));
D = [C*G; Xi];
delta = [f-C*b; zeros(size(Xi,1),1)];
Omega = blkdiag(Omega_f,eye(size(Xi,1)));
% Factor the imposed covariance BEFORE transmission through D's inverse.
% A weak driving shock can amplify one direction enormously. Thresholding
% eigenvalues of the resulting shock covariance would then erase unrelated,
% strictly positive soft variances.
L_Omega = blkdiag(psd_factor(Omega_f),eye(size(Xi,1)));
if isempty(D)
    U = zeros(0,0); V = eye(Nhor); singular = zeros(0,1); rank_D = 0;
else
    [U,S_D,V] = svd(D);
    q_D = min(size(D));
    singular = diag(S_D(1:q_D,1:q_D));
    rank_D = sum(singular > max(size(D))*eps(max(singular)));
end
Ur = U(:,1:rank_D); Vr = V(:,1:rank_D);
overdetermined = size(D,1) > Nhor;   % k > n*h: no exact solution, Eqs. (10)-(12) are the best approximate one (Sec. 2.2, App. B)
% Check attainability with an orthogonal projector, avoiding D*pinv(D)
% whose rounding error grows with the condition number of D.
if ~overdetermined && (norm(delta-Ur*(Ur'*delta)) > 1e-9*max(1,norm(delta)) || ...
        norm(L_Omega-Ur*(Ur'*L_Omega),'fro') > 1e-9*max(1,norm(L_Omega,'fro')))
    error('scenario:inconsistent','Scenario restrictions are not jointly attainable.');
end
sr = singular(1:rank_D);
mu_eps = Vr*((Ur'*delta)./sr);
L_restricted = Ur'*L_Omega;
L_eps = [Vr*bsxfun(@rdivide,L_restricted,sr), V(:,rank_D+1:end)];
if rank_D == Nhor && size(D,1) == Nhor
    % In a determined system, solve directly. This avoids the unnecessary
    % cancellation in V*diag(1./s)*U' for nearly triangular path constraints.
    mu_eps = D\delta;
    L_eps = D\L_Omega;
end
factor_target = [L_Omega,zeros(size(D,1),Nhor-rank_D)];
if any(~isfinite(mu_eps)) || any(~isfinite(L_eps(:))) || (~overdetermined && ( ...
        norm(D*mu_eps-delta)>1e-8*max(1,norm(delta)) || ...
        norm(D*L_eps-factor_target,'fro')>1e-8*max(1,norm(factor_target,'fro'))))
    error('scenario:numerical','Scenario constraints cannot be resolved to numerical tolerance.');
end
Sigma_eps = L_eps*L_eps';
mu_y = b+G*mu_eps;
L_y = G*L_eps;
Sigma_y = L_y*L_y';

%% 6. Scenario plausibility
% Compute the KL divergence of the scenario shocks from their unconditional distribution (Antolin-Diaz et al. 2021,
% equation (16)) and its calibration q (Section 3).
% In the orthogonal [Vr,Vnull] basis the covariance is block diagonal:
% diag(1./sr)*(L_restricted*L_restricted')*diag(1./sr), and I.
% Determine singularity before amplification and compute its log determinant
% from singular values of the factor, never eigenvalues of Sigma_eps.
restricted_sd = svd(L_restricted);
if numel(restricted_sd) < rank_D || any(restricted_sd == 0)
    det_Sigma_eps = 0; kl = Inf; q = 1;
else
    log_det = 2*(sum(log(restricted_sd))-sum(log(sr)));
    det_Sigma_eps = exp(log_det); % diagnostic only
    kl = max(0,0.5*(sum(L_eps(:).^2)+mu_eps'*mu_eps-Nhor-log_det));
    q = 0.5*(1+sqrt(-expm1(-2*kl/Nhor)));
end

%% 7. Conditional future paths
% Draw paths from the conditional Gaussian distribution (Antolin-Diaz et al. 2021, Algorithm 1).
y_draws = zeros(fhor,Nvar,0);
eps_draws = zeros(fhor,Nvar,0);
if Ndraws > 0
    z = randn(size(L_eps,2),Ndraws);
    eps_sim = repmat(mu_eps,1,Ndraws) + L_eps*z;
    y = repmat(mu_y,1,Ndraws)+L_y*z;
    eps_draws = permute(reshape(eps_sim,Nvar,fhor,Ndraws),[2 1 3]);
    y_draws = permute(reshape(y,Nvar,fhor,Ndraws),[2 1 3]);
end

%% 8. Output
% Collect moments, draws, restrictions and diagnostics in the result structure.
result.y_mean = reshape(mu_y,Nvar,fhor)';
result.y_mean_vector = mu_y;
result.y_cov = Sigma_y;
result.eps_mean = reshape(mu_eps,Nvar,fhor)';
result.eps_mean_vector = mu_eps;
result.eps_cov = Sigma_eps;
result.eps_factor = L_eps;
result.y_factor = L_y;
result.y_draws = y_draws;
result.eps_draws = eps_draws;
result.system = struct('b',b,'G',G);
result.constraints = struct('C',C,'f',f,'Omega_f',Omega_f, ...
    'Xi',Xi,'D',D,'delta',delta,'Omega',Omega);
result.shocks = struct('driving',driving,'non_driving',non_driving, ...
    'n_restrictions',size(Xi,1));
result.plausibility = struct('kl_divergence',kl,'q',q, ...
    'det_Sigma_eps',det_Sigma_eps,'dimension',Nhor, ...
    'direction','scenario_to_unconditional');
result.specification = scenario;
if isfield(scenario,'name'), result.name = scenario.name; else, result.name = ''; end
end

function L = psd_factor(A)
% Preserve every positive variance in a positive-definite input, regardless
% of scale. Only the semidefinite fallback removes numerical zero modes.
if isempty(A) || ~any(A(:))
    L = zeros(size(A,1),0); return
end
[L,flag] = chol(A,'lower');
if flag == 0, return; end
[V,E] = eig((A+A')/2); e = real(diag(E));
tol = size(A,1)*eps(max(abs(e)));
if any(e < -max(tol,1e-10*max(1,norm(A,2))))
    error('scenario:covariance','Imposed covariance must be positive semidefinite.');
end
active = e > tol;
L = V(:,active)*diag(sqrt(e(active)));
end
