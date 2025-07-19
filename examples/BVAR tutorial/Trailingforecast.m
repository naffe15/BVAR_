%%%%%%%%%%%
% GDP Nowcast
%%%%%%%%%%%

close all; clc; clear ;

addpath ../../cmintools/
addpath ../../bvartools/

% load the data
load DataMF
y = [GDP IPI HICP CORE Euribor1Y UNRATE];
lags = 6;

%% Estimate the MFVAR model
options.mf_varindex     = 1; 
options.K               = 1090;   
options.priors.name     = 'Minnesota';
options.noprint         = 1;
LastDataPoint           = find(T==2025)-1; % december 2024
bvarmf                  = bvar_(y(1:LastDataPoint,:),lags,options);

%% Construct the nowcast dataset
Tstart = 1;
YNowCast = nan(LastDataPoint  -Tstart+1 +3, size(y,2));
YNowCast(1:LastDataPoint,:) = y(Tstart:LastDataPoint,:);
yNowCast = repmat(YNowCast,1,1,3*5);
jj = 0; tmp = YNowCast;
for mm = 1 : 3 % month (jan, feb, mar)
    for vv = 2 : 6 % variable (IPI HICP CORE Euribor1Y UNRATE)
        jj = jj + 1;
        % creating a dataset with increasing amount of information
        yNowCast(:,:,jj) = tmp;
        % adding one more datapoint
        yNowCast(LastDataPoint -Tstart+1 + mm,vv,jj) = y(LastDataPoint+mm,vv);
        % storing the trailing data
        tmp = yNowCast(:,:,jj);
    end
end

%% Run the nowcast
NowOptions.fast_kf = 1;
NowOptions.noprint = 1;
out = nowcast_bvar(yNowCast, bvarmf, NowOptions);
% fast_kf=1 
% Elapsed time is 26.411855 seconds.
% fast_kf=0
% Elapsed time is 353.161598 seconds.

%% plot nowcast against actual

true = y(LastDataPoint + 3,1);
for dd = 1: size(yNowCast,3)
     GDP_now(:,dd) = squeeze(out.NowCast(LastDataPoint+3,:,dd));
end    
GDP_dn = prctile(GDP_now,5,1);
GDP_up = prctile(GDP_now,95,1);


tmp_str = 'mfvar_plt';
mkdir(tmp_str);
% string = IPI HICP CORE Euribor1Y UNRATE
figure,
x = 1:dd;
plot(dd,true,'kd','MarkerSize',12,'LineWidth',2);
hold on;
plot(x,median(GDP_now,1),'Color',[0.6 0.8 1],'LineWidth',1.5);
fill([x fliplr(x)], [GDP_up fliplr(GDP_dn)], [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
title('GDP - 2025Q1')
xlim([1 dd+0.1])
xlabel('data releases (#)');
ylabel('log levels');
grid on;
legend('actual','Median','90% CI');
ylim padded;
set(    gcf,'position' ,[50 50 900 650])
if strcmp(version('-release'),'2022b') == 0
    savefigure_pdf([tmp_str '\GDP_Now']);
end

figure,
x = 1:dd;
plot(dd,100*(true-y(LastDataPoint+3-12,1)),'kd','MarkerSize',12,'LineWidth',2);
hold on;
plot(x,100*(median(GDP_now,1)-y(LastDataPoint+3-12,1)),'Color',[0.6 0.8 1],'LineWidth',1.5);
fill([x fliplr(x)], [100*(GDP_up-y(LastDataPoint+3-12,1)) fliplr(100*(GDP_dn-y(LastDataPoint+3-12,1)))], [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
title('GDP annual growth rates - 2025Q1')
xlim([1 dd+0.1])
xlabel('data releases (#)');
ylabel('percent');
grid on;
legend('actual','Median','90% CI');
ylim padded;
set(    gcf,'position' ,[50 50 900 650])
if strcmp(version('-release'),'2022b') == 0
    savefigure_pdf([tmp_str '\GDP_Now_Growth']);
end

figure,
x = 1:dd;
plot(dd,100*(true-y(LastDataPoint+3-3,1)),'kd','MarkerSize',12,'LineWidth',2);
hold on;
plot(x,100*(median(GDP_now,1)-y(LastDataPoint+3-3,1)),'Color',[0.6 0.8 1],'LineWidth',1.5);
fill([x fliplr(x)], [100*(GDP_up-y(LastDataPoint+3-3,1)) fliplr(100*(GDP_dn-y(LastDataPoint+3-3,1)))], [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
title('GDP quarterly growth rates - 2025Q1')
xlim([1 dd+0.1])
xlabel('data releases (#)');
ylabel('percent');
grid on;
legend('actual','Median','90% CI');
ylim padded;
set(    gcf,'position' ,[50 50 900 650])
if strcmp(version('-release'),'2022b') == 0
    savefigure_pdf([tmp_str '\GDP_Now_GrowthQ']);
end


%%
return

%% Generate the first 5 datasets
 baseData = bvarmf.yinterpol(1:300, :);
% Specify the number of variables (columns) to use
max_vars = 5;
% Generation loop
for i = 1:max_vars
    % Empty matrix for nowcast variables
    nan_rows = NaN(3, size(baseData, 2)); 
    for j = 1:i
        %Fill the empty matrix
        nan_rows(1, j+1) = y(301, j+1);
    end
    % Save the dataset
    DataMFNC{i} = [baseData; nan_rows];
end

%% Generate the remaining datasets
% Matrix of indeces of values to add
to_add = [ ...
    302, 2;
    302, 3;
    302, 4;
    302, 5;
    302, 6;
    303, 2;
    303, 3;
    303, 4;
    303, 5;
    303, 6];
% Generation loop
for k = 1:size(to_add,1)
    r = to_add(k,1) - 300;
    c = to_add(k,2);
    % Fill the empty matrix
    nan_rows(r, c) = bvarmf.yinterpol(300 + r, c);
    % Save the dataset
    DataMFNC{5 + k} = [baseData; nan_rows];
end

%% Run the Kalman filter 
% Define the indeces of every dataset
index_list = {
    [2 2 0 0 0 0];  % DataMF1
    [2 2 2 0 0 0];  % DataMF2
    [2 2 2 2 0 0];  % DataMF3
    [2 2 2 2 2 0];  % DataMF4
    [2 2 2 2 2 2];  % DataMF5 
    [2 2 0 0 0 0];  % DataMF6
    [2 2 2 0 0 0];  % DataMF7
    [2 2 2 2 0 0];  % DataMF8
    [2 2 2 2 2 0];  % DataMF9
    [2 2 2 2 2 2];  % DataMF10
    [2 2 0 0 0 0];  % DataMF11
    [2 2 2 0 0 0];  % DataMF12
    [2 2 2 2 0 0];  % DataMF13
    [2 2 2 2 2 0];  % DataMF14
    [2 2 2 2 2 2];  % DataMF15
};
KFresults = cell(15, 100);
% Run the Kalman Filter for every dataset
% nowcast = nan(,length(DataMFNC),100);

for j = 1:15
    data = DataMFNC{j};
    KFoptions.index = [2
     0
     0
     0
     0
     0]; 

    % 1st variable is quarterly
    % KFoptions.mf_variables_kfindex = 1;   
    for i = 1:100%900:1000
        Phi   = bvarmf.Phi_draws(:,:,i);
        % Companion_matrix(1:bvarmf.N,:) = Phi(1:bvarmf.N*lags,:)';
        % test = (abs(eig(Companion_matrix)));
        % if any(test>1.0000000000001)
            KFoptions.initialCond = 1;
        % end
        Sigma = bvarmf.Sigma_draws(:,:,i);
        [logL, KFout] = kfilternan(Phi, Sigma, data, KFoptions);
        %KFresults{j, i-899} = KFout;
        KFresults{j, i} = KFout;
        % nowcast(:,i,j) = KFout.smoothSt_plus_ss(:,KFout.index_var(KFoptions.mf_variables_kfindex));
    end
end

%% Plot with estimation for every month
GDP_jan = zeros(15,100);
GDP_feb = zeros(15,100);
GDP_mar = zeros(15,100);
for j = 1:15
    for d = 1:100
        smoothed = KFresults{j,d}.smoothSt_plus_ss;
        GDP_jan(j,d) = smoothed(end-2,1);
        GDP_feb(j,d) = smoothed(end-1,1);
        GDP_mar(j,d) = smoothed(end,1);
    end
end
for j = 1:15
    % Sort the draws
    sorted_jan = sort(GDP_jan(j,:));
    sorted_feb = sort(GDP_feb(j,:));
    sorted_mar = sort(GDP_mar(j,:));
    % Relevant percentiles for each month
    low_jan(j)   = sorted_jan(5);
    high_jan(j)  = sorted_jan(95);
    med_jan(j)   = sorted_jan(50);
    low_feb(j)   = sorted_feb(5);
    high_feb(j)  = sorted_feb(95);
    med_feb(j)   = sorted_feb(50);
    low_mar(j)   = sorted_mar(5);
    high_mar(j)  = sorted_mar(95);
    med_mar(j)   = sorted_mar(50);
end
x = 1:15;
figure; hold on;
% January
fill([x fliplr(x)], [high_jan fliplr(low_jan)], [0.6 0.8 1], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(x, med_jan, 'b', 'LineWidth', 2);
% February
fill([x fliplr(x)], [high_feb fliplr(low_feb)], [0.8 1 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(x, med_feb, 'g', 'LineWidth', 2);
% March
fill([x fliplr(x)], [high_mar fliplr(low_mar)], [1 0.6 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(x, med_mar, 'r', 'LineWidth', 2);
xlabel('Updated datasets (#)');
ylabel('GDP Nowcast');
legend({'Jan CI','Jan median','Feb CI','Feb median','Mar CI','Mar median'}, 'Location','best');
title('GDP Nowcast Jan–Mar 2025');
grid on;


%% Plots for every month in one grid
x = 1:15;
figure;
% January
subplot(1,3,1);
fill([x fliplr(x)], [high_jan fliplr(low_jan)], [0.6 0.8 1], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(x, med_jan, 'b', 'LineWidth', 2);
title('GDP Nowcast – January 2025');
xlabel('Updated datasets (#)');
ylabel('GDP Nowcast');
grid on;
ylim padded;
% February
subplot(1,3,2);
fill([x fliplr(x)], [high_feb fliplr(low_feb)], [0.8 1 0.6], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(x, med_feb, 'g', 'LineWidth', 2);
title('GDP Nowcast – February 2025');
xlabel('Updated datasets (#)');
grid on;
ylim padded;
% March
subplot(1,3,3);
fill([x fliplr(x)], [high_mar fliplr(low_mar)], [1 0.7 0.7], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(x, med_mar, 'r', 'LineWidth', 2);
title('GDP Nowcast – March 2025');
xlabel('Updated datasets (#)');
grid on;
ylim padded;

%% Individual plots not in grid
x = 1:15;
% January
figure;
fill([x fliplr(x)], [high_jan fliplr(low_jan)], [0.6 0.8 1], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(x, med_jan, 'b', 'LineWidth', 2);
title('GDP Nowcast – January 2025');
xlabel('Updated datasets (#)');
ylabel('GDP Nowcast');
grid on;
legend('90% CI', 'Median');
ylim padded;
% February
figure;
fill([x fliplr(x)], [high_feb fliplr(low_feb)], [0.8 1 0.6], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(x, med_feb, 'g', 'LineWidth', 2);
title('GDP Nowcast – February 2025');
xlabel('Updated datasets (#)');
ylabel('GDP Nowcast');
grid on;
legend('90% CI', 'Median');
ylim padded;
% March
figure;
fill([x fliplr(x)], [high_mar fliplr(low_mar)], [1 0.7 0.7], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4); hold on;
plot(x, med_mar, 'r', 'LineWidth', 2);
title('GDP Nowcast – March 2025');
xlabel('Updated datasets (#)');
ylabel('GDP Nowcast');
grid on;
legend('90% CI', 'Median');
ylim padded;

%% Plot with real value comparison for March 2025
% Extract Real GDP value for March 2025
true_val = bvarmf.yinterpol(303, 1);
x = 1:15;
figure; hold on;
fill([x fliplr(x)], [high_mar fliplr(low_mar)], [1 0.7 0.7], ...
     'EdgeColor', 'none', 'FaceAlpha', 0.4);
plot(x, med_mar, 'r', 'LineWidth', 2);
% True value line
yline(true_val, '--k', 'LineWidth', 2);
xlabel('Updated datasets (#)');
ylabel('GDP Nowcast');
title('GDP Nowcast – March 2025');
legend({'90% CI', 'Median', 'True value'}, 'Location','best');
grid on;
ylim padded;
