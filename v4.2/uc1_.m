function [ytrend,ycycle] =uc1_(y,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates the UC model with correlated errors.

% Reference:
% Grant, A.L. and Chan, J.C.C. (2017). A Bayesian Model Comparison for
% Trend-Cycle Decompositions of Output, Journal of Money, Credit and Banking,
% 49(2-3): 525-552
% UC model
% y(t) = tau(t) + c(t), (1)
% where tau(t) is the trend and c(t) is the stationary, cyclical component.
% tau(t) = µ + tau(t-1) + utau(t)               (2)
% c(t)   = ph1 c(t-1) + phi2 c(t-2) + uc(t)     (3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = length(y);
time = 1:T;
% these priors are calibrated to the US real data.
% Prior mean µ
mu0         = .1;
% Prior variance µ
Vmu         = 1^2;
% Prior mean ?(0)
tau00       = 750;
% Prior variance tau(0)
Vtau0       = 100;
% Prior mean phi1 and phi2
phi0        = [0.5 0.2]';
invVphi     = speye(2);
% Prior scale uc(t)
sigc2_ub    = 0.5;
% Prior scale u?(t)
sigtau2_ub  = 1.5;
%
K       = 50000;
rho     = 0; %-.1;
rhoyes  = 0;% no correlation between trend and cycles 
ngrid   = 100; %500;

nsims   = 5000;    % MCMC parameter: number  of simulations
burnin  = 30000;  % MCMC parameter: number  of  burn-in

noprint     = 0; % print estimate
figg        = 1; %plot on
varname ={'Var1'};

if nargin>1
    if isfield(options,'figg') == 1
        figg = options.figg;
    end
    if isfield(options,'mu0') == 1
        % Prior mean µ
        mu0         = options.mu0;
    end
    if isfield(options,'Vmu') == 1
        % Prior variance µ
        Vmu         = options.Vmu;
    end
    if isfield(options,'tau00') == 1
        % Prior mean tau(0), ideally large value
        tau00       = options.tau00;
    end
    if isfield(options,'Vtau0') == 1
        % Prior variance tau(0)
        Vtau0       = options.Vtau0;
    end
    if isfield(options,'phi0') == 1
        % Prior mean phi1 and phi2
        phi0        = options.phi0;
        invVphi     = speye(2);
    end
    if isfield(options,'Vphi') == 1
        % Prior mean phi1 and phi2
        Vphi        = options.Vphi;
        invVphi     = inv(Vphi);
    end
    if isfield(options,'sigc2') == 1
        % Prior scale uc(t)
        sigc2_ub    = options.sigc2;
    end
    if isfield(options,'sigtau2') == 1
        % Prior scale u?(t)
        sigtau2_ub  = options.sigtau2;
    end
    if isfield(options,'rhoyes') == 1
        % correlation between cycles and trend
        rhoyes = options.rhoyes;
    end
    if isfield(options,'nsims') == 1
        nsims           = options.nsims;
    end
    if isfield(options,'burnin') == 1
        burnin          = options.burnin;
    end
    if isfield(options,'noprint')==1
        noprint = options.noprint;
    end
    if isfield(options,'time') ==1
        time = options.time;
        if length(time) ~= length(y)
            error('Time must be of the same size as y')
        end
    end
    if isfield(options,'varname') ==1
        varname = options.varname;
    end

    
end

%% Initialize the Markov chain
% Initial value for µ
mu      = mu0;%.1;
% Initial Values for phi1 and phi2
phi     = phi0;%[0.5 0.2]';
H       = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HH      = H'*H;
Hphi    = speye(T) - phi(1)*sparse(2:T,1:(T-1),ones(1,T-1),T,T) + ...
    - phi(2)*sparse(3:T,1:(T-2),ones(1,T-2),T,T);
% Initial value for tau(0)
tau0    = 1.0*y(1);
% Initial value for sigma_uc
sigc2   = sigc2_ub; 
% Initial value for sigma_utau
sigtau2 = sigtau2_ub;


%% prior evalution
pri_sigc2       = @(x) log(1/sigc2_ub);
pri_sigtau2     = @(x) log(1/sigtau2_ub);
pri_rho         = @(x) log(1/2);
count = 0;
tempphi = repmat(phi0',K,1) + (chol(invVphi\speye(2),'lower')*randn(2,K))';
for i=1:K
    phic = tempphi(i,:)';
    if sum(phic) < .99 && phic(2) - phic(1) < .99 && phic(2) > -.99
        count = count+1;
    end
end
phi_const = 1/(count/K);
prior = @(m,ph,sy,st,r,t0) -.5*log(2*pi*Vmu) -.5*(m-mu0)^2/Vmu ...
    -log(2*pi)+.5*log(det(invVphi))+log(phi_const)-.5*(ph-phi0)'*invVphi*(ph-phi0)...
    + pri_sigc2(sy) + pri_sigtau2(st) + pri_rho(r) ...
    -.5*log(2*pi*Vtau0) - .5*(t0-tau00)^2/Vtau0;

%disp('Starting MCMC for UCUR.... ');
start_time = clock;

%% initialize for storage
store_theta = zeros(nsims,7); % [mu, phi, sigc2, sigtau2, rho, tau0]
store_tau = zeros(nsims,T);
countphi = 0;

%rand('state', sum(100*clock) ); randn('state', sum(200*clock) );
rng('default' );

for isim = 1:nsims+burnin
    
    %% sample tau
    alp = H\(mu+[tau0;sparse(T-1,1)]);
    a = -rho*sqrt(sigc2/sigtau2)*(H*alp);
    B = Hphi + rho*sqrt(sigc2/sigtau2)*H;
    tmpc = 1/((1-rho^2)*sigc2);
    Ktau = HH/sigtau2 + tmpc*(B'*B);
    %    size(Ktau)
    %    size(HH)
    %    size(alp)
    %    size(sigtau2)
    %    size(tmpc)
    %    size(B)
    %    size(Hphi)
    %    size(y)
    tauhat = Ktau\(HH*alp/sigtau2 + tmpc*B'*(Hphi*y-a));
    tau = tauhat + chol(Ktau,'lower')'\randn(T,1);
    
    
    
    %% sample phi
    e = y-tau;
    Xphi = [[0;e(1:T-1)] [0;0;e(1:T-2)]];
    tmpc = 1/((1-rho^2)*sigc2);
    Kphi = invVphi + tmpc*(Xphi'*Xphi);
    phihat = Kphi\(invVphi*phi0 ...
        + tmpc*Xphi'*(e-rho*sqrt(sigc2/sigtau2)*H*(tau-alp)));
    flag = 0; count = 0;
    while flag == 0 && count < 100
        phic = phihat + chol(Kphi,'lower')'\randn(2,1);
        if sum(phic) < .99 && phic(2) - phic(1) < .99 && phic(2) > -.99
            phi = phic;
            flag = 1;
            countphi = countphi + 1;
        end
        count = count + 1;
    end
    Hphi = speye(T) - phi(1)*sparse(2:T,1:(T-1),ones(1,T-1),T,T) + ...
        - phi(2)*sparse(3:T,1:(T-2),ones(1,T-2),T,T);
    
    %% sample sigc2
    u = [Hphi*(y-tau) (tau-[tau0; tau(1:end-1)])-mu];
    c1 = sum(u(:,1).^2);
    c2 = u(:,1)'*u(:,2);
    c3 = sum(u(:,2).^2);
    gy = @(x) -T/2*log(x) - 1./(2*(1-rho^2)*x).*(c1-2*rho*sqrt(x/sigtau2)*c2...
        + rho^2*x/sigtau2*c3);
    sigc2grid = linspace(rand/100,sigc2_ub-rand/100,ngrid);
    logpsigc2 = gy(sigc2grid) + pri_sigc2(sigc2grid);
    psigc2 = exp(logpsigc2-max(logpsigc2));
    psigc2 = psigc2/sum(psigc2);
    cumsumy = cumsum(psigc2);
    sigc2 = sigc2grid(find(rand<cumsumy, 1 ));
    
    %% sample sigtau2
    gtau = @(x) -T/2*log(x) - c3./(2*x) ...
        - 1/(2*(1-rho^2)*sigc2)*(c1-2*rho*sqrt(sigc2./x)*c2 + rho^2*sigc2./x*c3);
    sigtau2grid = linspace(rand/100,sigtau2_ub-rand/100,ngrid);
    logpsigtau2 = gtau(sigtau2grid) + pri_sigtau2(sigtau2grid);
    psigtau2 = exp(logpsigtau2-max(logpsigtau2));
    psigtau2 = psigtau2/sum(psigtau2);
    cumsumtau = cumsum(psigtau2);
    sigtau2 = sigtau2grid(find(rand<cumsumtau, 1 ));
    
    %% sample rho
    
    if rhoyes==1
        grho = @(x) -T/2*log(1-x.^2) + ...
            - 1./(2*sigc2*(1-x.^2)).*(c1-2*x*sqrt(sigc2/sigtau2)*c2 ...
            + x.^2*sigc2/sigtau2*c3);
        rhogrid = linspace(-.98-rand/100,.98+rand/100,ngrid);
        logprho = pri_rho(rhogrid) + grho(rhogrid);
        prho = exp(logprho-max(logprho));
        prho = prho/sum(prho);
        cumsumrho = cumsum(prho);
        rho = rhogrid(find(rand<cumsumrho, 1));
    else
        rho=0.0;
    end
    %% sample mu and tau0
    uc = Hphi*(y-tau);
    Xdel = [ones(T,1) H\ones(T,1)];
    Kdel = diag([1/Vtau0  1/Vmu]) + Xdel'*HH*Xdel/((1-rho^2)*sigtau2);
    delhat = Kdel\([tau00/Vtau0; mu0/Vmu] + 1/((1-rho^2)*sigtau2)*Xdel'*HH*...
        (tau - rho*sqrt(sigtau2/sigc2)*(H\uc)));
    del = delhat + chol(Kdel,'lower')'\randn(2,1);
    tau0 = del(1);
    mu = del(2);
    
    %   if ( mod( isim, 10000 ) ==0 )
    %       disp(  [ num2str( isim ) ' loops... ' ] )
    %   end
    
    if isim > burnin
        isave = isim - burnin;
        store_tau(isave,:) = tau';
        store_theta(isave,:) = [mu phi' sigc2 sigtau2 rho tau0];
    end
end


disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );

% if cp_ml
%     start_time = clock;
%     disp('Computing the marginal likelihood.... ');
%     [ml,mlstd] = ml_UCUR(y,store_theta,prior,M);
%     disp( ['ML computation takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
% end

tauhat = median(store_tau)';
tauCI = quantile(store_tau,[.1 .9])';
thetahat = mean(store_theta)';
thetastd = std(store_theta)';

% print  estimates
if  noprint==0
    fprintf('\n');
    fprintf('Parameter   | Posterior mean (Posterior std. dev.):\n');
    fprintf('mu          | %.2f (%.2f)\n', thetahat(1), thetastd(1));
    fprintf('phi_1       | %.2f (%.2f)\n', thetahat(2), thetastd(2));
    fprintf('phi_2       | %.2f (%.2f)\n', thetahat(3), thetastd(3));
    fprintf('sigma^2_c   | %.2f (%.2f)\n', thetahat(4), thetastd(4));
    fprintf('sigma^2_tau | %.2f (%.2f)\n', thetahat(5), thetastd(5));
    fprintf('rho         | %.2f (%.2f)\n', thetahat(6), thetastd(6));
end

ycycle=y-tauhat;
ytrend=tauhat;

% % print  marginal  likelihood
% if cp_ml
%     fprintf('\n');
%     fprintf('log marginal likelihood: %.1f (%.2f)\n', ml, mlstd);
% end

%% plot estimated components
if figg==1
    %figure;
    subplot(2,1,1);
    hold on
    title(varname)
    %   plot(tid,tauhat, 'LineWidth',1,'Color','blue');
    plot(time,y,'r-','LineWidth',2);
    plot(time,ytrend, 'k--','LineWidth',2); axis tight;
    hold off
    box off; %xlim([tid(1)-1 tid(end)+1]);
    legend('Level','Trend','Location','NorthWest');
    subplot(2,1,2);
    tmpy = repmat(y,1,2)-tauCI;
    hold on
    %    plotCI(tid,tmpy(:,1),tmpy(:,2));
    %    plot(tid, (y-tauhat), 'LineWidth',1,'Color','blue');
    %    plot(tid, zeros(T,1),'-k','LineWidth',1);
    plot(time,ycycle, 'b','LineWidth',2); axis  tight;
    %    plot(tid, zeros(T,1),'-k','LineWidth',1);
    hold off
    box off; %xlim([tid(1)-1 tid(end)+1]); ylim([-10 8])
    legend('Cycle','Location','NorthWest');
    set(gcf,'Position',[100 100 800 300])  %% 1 by 2 graphs
    pause
    close all;
end



