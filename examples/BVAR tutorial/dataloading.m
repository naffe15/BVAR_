%% Forecast Databse
clear all
%EA: monthly database with various macro variables
datafile = 'data_macro.xlsx';
sheet    = 'logs';

[num,txt,raw] = xlsread(datafile,sheet);

% generate a column vector with the T span
timev  = num(1:end,1:2);
% this is MONTHLY data
% convention 2012m1 = 2012.00 and 2012m12 = 2012+11/12

T   = timev(1,1) + (timev(1,2)-1)/12 : 1/12 : timev(end,1) + (timev(end,2)-1)/12;

% remove the first two comumns of numbers (T span)
num(:,1:2) = [];

% generate a cell vector with the varables names
varnames  = txt(1,3:end);


save Data T num varnames
for jj = 1 : length(varnames)
    eval([varnames{jj} '=num(:,' num2str(jj) ');'])
    eval([ 'save Data ' varnames{jj} ' -append'])
end

%% MP shocks Database
clear all
%US: monthly GK dataset
datafile = 'VAR_data';
sheet    = 'VAR_data';

[num,txt,raw] = xlsread(datafile,sheet);

% generate a column vector with the T span
timev  = num(1:end,1:2);
% this is MONTHLY data
% convention 2012m1 = 2012.00 and 2012m12 = 2012+11/12

T   = timev(1,1) + (timev(1,2)-1)/12 : 1/12 : timev(end,1) + (timev(end,2)-1)/12;

% remove the first two comumns of numbers (T span)
num(:,1:2) = [];

% generate a cell vector with the varables names
varnames  = txt(1,3:end);

save DataGK T num varnames
for jj = 1 : length(varnames)
    eval([varnames{jj} '=num(:,' num2str(jj) ');'])
    eval([ 'save DataGK ' varnames{jj} ' -append'])
end

%% Mixed Freq Database
clear all

ploton =1;
% monthly database with various macro variables
[num,txt,~] = xlsread('data_macro.xlsx','mixed_freq');

% the standard interpolation only uses the individual TS info while the VAR
% more info. this is bad for end of sample nowcasts. to see this try to
% remove the last data points as below
% [num,txt,~] = xlsread('data_macro.xlsx','mixed_freq','A1:L186');

% generate a column vector with the T span
timev  = num(1:end,1:2);
% this is MONTHLY data change for Q data
T   = timev(1,1) + timev(1,2)/12 : 1/12 : timev(end,1) + timev(end,2)/12;
%
% remove the first two comumns of numbers (T span)
num(:,1:2) = [];

% generate a cell vector with the varables names
varnames  = txt(1,3:end);

ordering = {'GDP','IPI','HICP','CORE','Euribor1Y','UNRATE'};
% ordering = {'GDP','IPI'};
% this function extract the index of ordering in the original database and
% it generates index_XXX for the ordering in the VAR
[sindex] = OrderingIndexes(ordering,varnames);

% select the dataset with the prespecifed variable selection in ordering
y = num(:,sindex);

if ploton ==1
    figure('Name','Data')
    set(gcf,'position' ,[50 50 800 650])
    for jj=1:length(ordering)
        subplot(3,3,jj)
        plot(T,y(:,jj),'Linewidth',2,'marker','+')
        title(ordering{jj})
        axis tight
    end
end

save DataMF T num varnames
for jj = 1 : length(varnames)
    eval([varnames{jj} '=num(:,' num2str(jj) ');'])
    eval([ 'save DataMF ' varnames{jj} ' -append'])
end


%% Rolling Window Quarterly Database

clear all
datafile = 'Qdata.xls';
sheet    = 'Sheet1';

[num,txt,raw] = xlsread(datafile,sheet);

% generate a column vector with the T span
timev  = datevec(txt(2:end,1),'mm/dd/yyyy');
% this is Q data
% convention 2012q1 = 2012.00 and 2012q4 = 2012+3/4

T   = timev(1,1) + (timev(1,2)-1)/12 : 1/4 : timev(end,1) + (timev(end,2)-1)/12;

% generate a cell vector with the varables names
varnames  = txt(1,2:end);


save DataQ T num varnames
for jj = 1 : length(varnames)
    eval([varnames{jj} '=num(:,' num2str(jj) ');'])
    eval([ 'save DataQ ' varnames{jj} ' -append'])
end

%% FAVAR Database (Stock_Watson type database)

clear all
datafile = 'Stock_Watson_GR.xls';
sheet    = 'FAVAR';

[num,txt,raw] = xlsread(datafile,sheet);

% generate a column vector with the T span
timev  = datevec(txt(2:end,1),'mm/dd/yyyy');
% this is Q data
% convention 2012q1 = 2012.00 and 2012q4 = 2012+3/4

T   = timev(1,1) + (timev(1,2)-1)/12 : 1/4 : timev(end,1) + (timev(end,2)-1)/12;

% generate a cell vector with the varables names
varnames    = txt(1,2:end);
% fast moving variable (interest rate)
varnames_y1 = {'TBILL3M'};
% slow moving variables
% variables in log levels 
varnames_y2 = setdiff(varnames,varnames_y1);

y1      = num(:,strmatch(varnames_y1,varnames));
[~,gg]  = ismember(varnames_y2,varnames);
yy2      = num(:,gg);

ff=1;
y1(1:ff) = []; 
T(1:ff)  = [];
% construct YoY growth rate 
y2      = 100*(yy2(1+ff:end,:) - yy2(1:end-ff,:));

save DataFAVAR T y1 y2 varnames_y1 varnames_y2

%% EA Unconventional MP shocks

clear all 
close all
clc

% EA MP events database
% Altavilla et al. / Journal of Monetary Economics 108 (2019) 162–179
datafile = 'Dataset_EA-MPD.xlsx';
sheet    = 'Monetary Event Window';
[num,txt,raw] = xlsread(datafile,sheet);

% generate a cell vector with the varables names
varnames  = txt(1,2:end);
ordering = {'DE10Y','DE2Y'};
[sindex] = OrderingIndexes(ordering,varnames);

% generate a column vector with the T span
timev  = datevec(txt(2:end,1),'mm/dd/yyyy');
T      = 1999 : 1/12 : 2019+11/12;
% select the 10y and 2y German Govt Bonds Yields (longer time span)
DE10Y  = zeros(length(T),1);
DE2Y   = zeros(length(T),1);


tt = 0;
for yy = 1999 : 2019
    for mm = 1:12    
        tt = 1 + tt;
        index = find(sum([yy mm]==timev(:,1:2),2)==2);
        if isempty(index) == 0
            % if more than one MP event in the same month sum the surprises
            DE10Y(tt,1) = sum(num(index,sindex(1)));
            DE2Y(tt,1)  = sum(num(index,sindex(2)));            
        end
    end
end

datafile = 'dataEAm';
sheet    = 'dataEAm';
[num0,~,~] = xlsread(datafile,sheet);
HICP   = num0(:,1); CORE   = num0(:,2);	IPI    = num0(:,3); 
LTR10Y = num0(:,4); UNR = num0(:,5); Euribor1Y = num0(:,6); %PCOMM = num0(:,7);

save DataEAMP T HICP IPI LTR10Y DE10Y DE2Y CORE UNR Euribor1Y

