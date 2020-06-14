function [output] = ols_reg(Y,X,options)

% this function computes the ols coefficients of Y on X
% Y TxN
% X Txk (interepct at the end)
% Y = XB + E

nindex      = find(sum(isnan([Y,X]),2)==0);
index       = find(sum(isnan([Y,X]),2)>0);
Y(index,:)  = [];
X(index,:)  = [];
[N,K]       = size(X);
robust_se_  = 0;
L           = round(N/4);

if nargin > 2
    if isfield(options,'robust_se_') == 1
        robust_se_ = options.robust_se_;
    end
    if isfield(options,'L') == 1
        L = options.L;
    end
end


ny = size(Y,2);

XX      = (X'*X);
iXX     = XX\eye(K);
Bols    = X\Y;                  % Bols   = (X'*X)\(X'*Y);
err     = Y - X*Bols;
se      = nan(K,ny);
Sols    = zeros(K*ny);

if robust_se_==2 % NW Robust SE
    %-----------------------------------------------------------------------%
    Serror     = (err'*err);
    nwWeights  = (L+1-(1:L))./(L+1);
    for j= 1 : L
        G      = (err(j+1 : N, :)'*err(1 : N-j , :));
        Serror = Serror + nwWeights(j) * (G + G');
    end
    Serror = Serror/(N-K);
    Sols   = kron(diag(diag(Serror)),iXX);
    se     = reshape(sqrt(diag(Sols)), K, ny);
        
elseif robust_se_ == 1 % Hamilton (1994), Ch 10 pag 282, eq (10.5.20)
    %-----------------------------------------------------------------------%
    for vv = 1: ny % equation by equation
        u= err(:,vv); 
        errs=X.*u;
        V0 = [errs'*errs] / N ; %regular weighting matrix
        for ind_i = (1:L)
            S = errs(1:N-ind_i,:)'*errs(1+ind_i:N,:) / N;
            V0 = V0 + (1 - ind_i/(L+1))*(S + S');
        end
        %     D    =   inv((X'*X)/N);
        %     varb = 1/N*D*V1*D;
        Solsj =  N * iXX * V0 * iXX;
        Sols((vv-1)*K+1 : vv*K, (vv-1)*K+1 : K*vv ) = Solsj;
    end
    Serror = (err'*err)/(N-K); % not sure this is the correct Coavariance of the shocks. 
    se     = reshape(sqrt(diag(Sols)), K, ny);
    
    
elseif robust_se_ == 5 % Matlab HAC function 
    %-----------------------------------------------------------------------%
    if  exist('hac') ==2
        % find constant
        index = find(sum(diff(X,2),1)~=0);
        indexC = find(sum(diff(X,2),1)==0);
        for vv = 1: ny
%             [EstCoeffCov0,se0,~] = hac(X(:,index),Y(:,ny),'display','off','type','HC'); 
            [EstCoeffCov0,se0,~] = hac(X(:,index),Y(:,ny),'display','off','bandwidth','AR1OLS'); % remove the intercept
            se(index,vv)  = se0(2:end);
            se(indexC,vv) = se0(1); %intercept
            Sols((vv-1)*K+1 : vv*K, (vv-1)*K+1 : K*vv ) = EstCoeffCov0; % order is not correct
        end
        Serror  = 1/(N-K)*(err'*err);
%         output.TtestRobust  = coeff./se;
%         output.pvalueRobust = tpdf(output.TtestRobust,N-K);
    else
        error('Matlab Econ Toolbox missing')
    end
else
    Serror  = 1/(N-K)*(err'*err);
    Sols    = kron(diag(diag(Serror)),iXX);
    se      = reshape(sqrt(diag(Sols)), K, ny);
end

Ttest  = Bols./ reshape(sqrt(diag(Sols)), K, ny);
ESS    = diag((X*Bols - mean(Y))'*(X*Bols - mean(Y)));
RSS    = diag(err' * err);
TSS    = diag((Y - mean(Y))' * (Y - mean(Y)));
R2     = ones(length(ESS),1) - RSS ./ TSS;
adjR2  = ones(length(ESS),1) - (ones(length(ESS),1) - R2)*(N-1)/(N-K);
Ftest  = ESS/(K-1) ./ diag(Serror);

% for v = 1 : ny
% %     output.logl(v,1)   = -N/2*log(2*pi*Serror(v,v)) - RSS(v,1)/(2*Serror(v,v));
%     [output.AIC(v,1), output.SIC(v,1), output.HQIC(v,1)] = IC(output.logl(v,1), N, K);
% end

output.beta   = Bols;                   % OLS estimator
output.error  = err;                   % (TxN) matrix of Residuals
output.e_ols  = err;                   % (TxN) matrix of Residuals
output.Serror = Serror;                 % Covariance matrix of Residuals
output.Sigma_ols = Serror;                 % Covariance matrix of Residuals
output.Sols   = Sols;                   % Covariance matrix of Bols
output.Ttest  = Ttest;                  % t-statistics
output.pvalue = tpdf(Ttest,N-K);        % p-value
output.Ftest  = Ftest;                  % F-test
output.R2     = R2;                     % R2
output.adjR2  = adjR2;                  % Adjusted R2
output.yfit   = X*Bols;                 % Fitted Values
output.N      = N;                      % # of observation used
output.K      = K;                      % # of regressors
output.nindex = nindex;                 % index of missing observation
output.index  = index;                  % index of observation used
output.se     = se;
output.XX     = XX;
output.X      = X;
output.Y      = Y;

for v = 1 : ny
%     output.logl(v,1)   = -N/2*log(2*pi*Serror(v,v)) - RSS(v,1)/(2*Serror(v,v));
    [output.AIC(v,1), output.SIC(v,1), output.HQIC(v,1)] = IC(output, N, K);
end

% if nargin > 2
%     [EstCov,se,coeff] = hac(X(:,1:end-1),Y,'display','off'); % remove the intercept
%     output.serob        = se(2:end);
%     output.serob(end+1) = se(1); %intercept
%     output.TtestRobust  = coeff./se;
%     output.pvalueRobust = tpdf(output.TtestRobust,N-K);
% end


