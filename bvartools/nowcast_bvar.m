function out = nowcast_bvar(DataMFNC, BVAR, options)

fast_kf             = 0;
KFoptions           = BVAR.KFoptions;
NumbDrawsFromPost   = 100;
lags                = BVAR.lags;
fhor                = BVAR.fhor;
forecast_data       = BVAR.forecast_data;

[T, ny, nD]         = size(DataMFNC);

if nargin>2
    % TBA
    if isfield(options,'NumbDrawsFromPost')==1
        NumbDrawsFromPost =  options.NumbDrawsFromPost;
    end
    if isfield(options,'fast_kf')==1
        fast_kf = options.fast_kf;
    end
    if isfield(options,'noprint')==1
        KFoptions.noprint = options.noprint;
    end
end
index = randi(BVAR.ndraws,NumbDrawsFromPost,1);
NowCast = nan(T, ny, nD);
Companion_matrix = diag(ones(ny*(lags-1),1),-ny);
frcsts_no_shock = nan(fhor, ny, NumbDrawsFromPost); 
frcsts_with_shocks = frcsts_no_shock;

if fast_kf
    % find the last full raw wihtout nans
    tmp = find(sum(isnan(DataMFNC(:,:,1)),2)==0);
    To = tmp(end);
    y = DataMFNC(1:To,:,1);
     for i = 1:NumbDrawsFromPost
        Phi   = BVAR.Phi_draws(:,:,index(i));
        Sigma = BVAR.Sigma_draws(:,:,index(i));
        [ ~, KFzero] = kfilternan(Phi, Sigma, y, KFoptions);
        pZero(:,:,i) = KFzero.CovSt(:,:,end);
        aZero(:,i)   = KFzero.smoothSt(end,:)';
        forecast_data.initval = KFzero.smoothSt_plus_ss(end-lags+1:end,KFzero.index_var);%yfill(end-lags+1:end, :, d);
        [frcsts_no_shock(:,:,i),frcsts_with_shocks(:,:,i)] = forecasts(forecast_data,Phi,Sigma,fhor,lags);
     end

end

wb = waitbar(0, 'Computing the nowcast with the MF BVAR');
KK = NumbDrawsFromPost*nD;
dd = 0;
for j =  1: nD
    data = DataMFNC(:,:,j);
    for i = 1:NumbDrawsFromPost
        dd= dd + 1;
        Phi   = BVAR.Phi_draws(:,:,index(i));
        Sigma = BVAR.Sigma_draws(:,:,index(i));
        Companion_matrix(1:ny,:) = Phi(1:ny*lags,:)';
        test = (abs(eig(Companion_matrix)));
        if any(test>1.0000000000001)
            KFoptions.initialCond = 1;
        end
        if fast_kf==1
            KFoptions.initialCond = 2;
            KFoptions.pZero = pZero(:,:,i);
            KFoptions.aZero = aZero(:,i);
            [~, KFout] = kfilternan(Phi, Sigma, data(To+1:end,:), KFoptions);
            NowCast(1:To,i,j) = KFzero.smoothSt_plus_ss(1:To,KFout.index_var(KFoptions.mf_varindex));
            NowCast(1+To:end,i,j) = KFout.smoothSt_plus_ss(:,KFout.index_var(KFoptions.mf_varindex));
        else
            [~, KFout] = kfilternan(Phi, Sigma, data, KFoptions);
            NowCast(:,i,j) = KFout.smoothSt_plus_ss(:,KFout.index_var(KFoptions.mf_varindex));
        end
        waitbar(dd/KK, wb);
    end
end
close(wb);

out.NowCast = NowCast;
out.options = options; 
out.ForeCast.no_shock = [out.NowCast; squeeze(frcsts_no_shock(:,KFoptions.mf_varindex,:))] ;
out.ForeCast.with_shocks = [out.NowCast; squeeze(frcsts_with_shocks(:,KFoptions.mf_varindex,:))] ;