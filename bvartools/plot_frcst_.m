function plot_frcst_(frcsts,y,time,options)
%PLOT_FRCST_ Fan chart of the forecasts, one subplot per variable.
%
%   PLOT_FRCST_(frcsts,y,time)
%   PLOT_FRCST_(frcsts,y,time,options)
%
%   Plots the data followed by the median forecast and credible bands. See
%   the forecasting section of the Hitchhiker's Guide.
%
%   Inputs
%     frcsts   fhor x ny x K   forecast draws, e.g. BVAR.forecasts.no_shocks
%     y        T x ny          data
%     time     T x 1           dates of the data
%     options  struct (optional)
%       .varnames        cell, one name per variable
%       .time_start      first date plotted (must be in time)
%       .conf_sig        coverage of the band (default 0.68)
%       .conf_sig_2      coverage of a second band
%       .nplots          1 x 2   subplot grid
%       .add_frcst       (T+fhor) x ny   additional data and forecasts
%       .order_transform 1 x ny  transformation per variable: 0 none, 1
%                        period change, 4 or 12 the 100 x 4- or 12-period
%                        change, 100 the 100 x period change
%       .saveas_dir, .saveas_strng   folder and name; saved only if both
%                                    are given
%
%%gda
%   Structural scenarios: frcsts can also be the output of BVAR_ run with
%   options.scenario, or a cell array with one output per scenario. One
%   forecast figure per scenario (blue credible sets, target paths in red)
%   and one figure with the plausibility q of the scenarios. The name of each
%   scenario is options.scenario.name of its run of BVAR_ (default 'Scenario i').
%   Options: the options above and .add_frcst_label (legend of .add_frcst) and
%   .fig_name. Defaults: .conf_sig 0.40, .conf_sig_2 0.68,
%   .time_start 24 periods before the end of the sample, .nplots with 3
%   columns, .add_frcst_label 'Unconditional (median)'.
%%gda
%   See also BVAR_, FORECASTS, SHADE.

%   Filippo Ferroni, 6/1/2015.

%%gda
% structural scenarios: one forecast figure per scenario and a figure with the plausibility q
if isstruct(frcsts) || iscell(frcsts)
    if isstruct(frcsts)
        frcsts = {frcsts};
    end
    if nargin < 4
        options = struct();
    end
    nscn  = numel(frcsts);
    nvr   = size(y,2);
    names = arrayfun(@(i) sprintf('Scenario %d',i),1:nscn,'UniformOutput',false);
    for i = 1 : nscn
        if isfield(frcsts{i}.forecasts.scenario,'name')
            names{i} = frcsts{i}.forecasts.scenario.name;
        end
    end
    if ~isfield(options,'nplots'),          options.nplots          = [ceil(nvr/3) min(nvr,3)];       end
    if ~isfield(options,'time_start'),      options.time_start      = time(max(1,length(time)-23));   end
    if ~isfield(options,'conf_sig'),        options.conf_sig        = 0.40;                           end
    if ~isfield(options,'conf_sig_2'),      options.conf_sig_2      = 0.68;                           end
    if ~isfield(options,'add_frcst_label'), options.add_frcst_label = 'Unconditional (median)';       end
    for i = 1 : nscn
        opt          = options;
        opt.scenario = frcsts{i}.forecasts.scenario;
        opt.fig_name = ['Structural scenario: ' names{i}];
        plot_frcst_(frcsts{i}.forecasts.conditional,y,time,opt);
    end
    % plausibility: q goes from 0.5 (same as the unconditional forecast) to 1 (implausible);
    % a hard scenario has q = 1 in every draw
    figure('name','Structural scenario: plausibility (KL divergence and q)','position',[80 80 1400 700]);
    for i = 1 : nscn
        q = frcsts{i}.forecasts.plausibility.q;
        subplot(ceil(nscn/3),min(nscn,3),i); hold on
        histogram(q,'BinEdges',linspace(0.5,1,21),'Normalization','probability');
        title(sprintf('%s\nq: med %.3f, [%.3f, %.3f]',names{i},median(q),min(q),max(q)), ...
            'Interpreter','none','FontSize',9);
        xlabel('q (always in [0.5, 1])'); ylabel('fraction of draws'); xlim([0.5 1]);
    end
    return
end
%%gda
if length(time) ~= length(y)
    error('Mismatch between time and in sample data (y)');
end
if size(time,2)>size(time,1)
    time=time';
end

hor     = size(frcsts,1);
nvar    = size(frcsts,2);
ndraws  = size(frcsts,3);
%%gda
% scenario: drop the draws with non-finite values (sort would move them to the top and shift the quantiles)
if nargin == 4 && isfield(options,'scenario')
    valid = reshape(all(all(isfinite(frcsts),1),2),[],1);
    frcsts = frcsts(:,:,valid);
    ndraws = size(frcsts,3);
end
scenario_yes = 0;
fig_name     = 'Forecasts';
add_label    = 'Additional series';
%%gda
nplots  = [ceil(sqrt(nvar)) ceil(sqrt(nvar))];
savefig_yes = 0;
conf_sig    = 0.68;
add_frcst_yes = 0;
add_multiple_bands_yes = 0;
fnam_dir    = '.';
fnam_suffix = 'frcsts';
trasf_yes   =  0;

% retreive the frequency
integerTest = ~mod(time,1);
indx        = find(integerTest==1);
frq         = indx(2)-indx(1);

if frq == 12 % monthly
       timefor = time(end) + 1/12 : 1/12 : time(end) + hor/12;
       time = [time; timefor'];
elseif frq == 4 % quarterly
       timefor = time(end) + 1/4 : 1/4 : time(end) + hor/4;
       time = [time; timefor'];
elseif frq == 1 % annual
       timefor = time(end) + 1 : 1 : time(end) + hor;
       time = [time; timefor'];
elseif frq == 48 % weekly
       timefor = time(end) + 1/48 : 1/48 : time(end) + hor/48;
       time = [time; timefor'];
else
    error('Frequency not defined.')
end
% 
% if strmatch(freq,'m') == 1 %#ok<*MATCH2>
%        timefor = time(end) + 1/12 : 1/12 : time(end) + hor/12;
%        time = [time; timefor'];
% elseif strmatch(freq,'q') == 1
%        timefor = time(end) + 1/4 : 1/4 : time(end) + hor/4;
%        time = [time; timefor'];
% elseif strmatch(freq,'a') == 1
%        timefor = time(end) + 1 : 1 : time(end) + hor;
%        time = [time; timefor'];
% else
%     error('You need to provide a frequency: ''m'', ''q'' or ''a''.')
% end
time_start = 1;
time_end = length(time);

if nargin < 4
    disp('You did not provided names for variables.')
    disp('I call them Var 1, Var 2, ... ')
    for v = 1 : nvar
        eval(['varnames{'   num2str(v) '} =  ''Var  ' num2str(v) ''';'])
    end
else    
    if isfield(options,'time_start') ==1 
        time_start = find(options.time_start==time);
        if isempty(time_start) ==1 
            error('''time_start'' is not included in ''T''.')
        end
    end
    if isfield(options,'order_transform') ==1 
        trasf_yes       = 1;
        order_trasform  = options.order_transform;
        if length(order_trasform) ~= nvar
            error('Mismatch between the ''order_transform'' size and # of variables.')
        end
    end
    if isfield(options,'varnames') ==1
        varnames = options.varnames;
        if length(varnames) ~= nvar
            error('Mismatch between the # varnames and # of variables to plot')
        end
    else
        disp('You did not provided names for the endogenous variables.')
        disp('I call them Var 1, Var 2, ...')
        for v = 1 : nvar
            eval(['varnames{'   num2str(v) '} =  ''Var  ' num2str(v) ''';'])
        end
    end
    if isfield(options,'nplots') ==1 
        nplots = options.nplots;
    end
    if isfield(options,'saveas_strng') ==1 
        savefig_yes = 1;
        % setting the names of the figure to save
        fnam_suffix = [ fnam_suffix options.saveas_strng ];
    end
    if isfield(options,'saveas_dir') == 1
        savefig_yes = 1;
        % setting the folder where to save the figure
        fnam_dir = options.saveas_dir;
        if exist(fnam_dir,'dir') == 0
            mkdir(fnam_dir)
        end
    end
    if isfield(options,'add_frcst') ==1
        add_frcst = options.add_frcst;
        if size(add_frcst) ~= [length(time), nvar]
            error('The ''add_frcst'' dimensions must be in-sample + output-of-sample length and the # of variables to plot')
        end
        add_frcst_yes = 1;
    end    
    if isfield(options,'conf_sig') ==1
        conf_sig = options.conf_sig;
    end
    if isfield(options,'conf_sig_2') ==1
        if options.conf_sig_2 < conf_sig
            error('Additional confidence bands should be larger than ''options.conf_sig''.')
        end
        add_multiple_bands_yes = 1;
        sort_idx_2   = round((0.5 + [-options.conf_sig_2, options.conf_sig_2, 0]/2) * ndraws);
    end
%%gda
    if isfield(options,'add_frcst_label'), add_label = options.add_frcst_label; end
    if isfield(options,'fig_name'),        fig_name  = options.fig_name;        end
    if isfield(options,'scenario')
        % BVAR.forecasts.scenario: a row of C with a single nonzero in column (h-1)*nvar+v
        % pins variable v at forecast period h
        scenario_yes = 1;
        C   = options.scenario.C;
        f   = options.scenario.f(:);
        one = sum(C~=0,2) == 1;
        if ~all(one)
            warning('%d scenario restrictions combine several variables and are not plotted.',sum(~one));
        end
        [~,col]  = max(abs(C(one,:)),[],2);
        tgt      = nan(nvar,hor);
        tgt(col) = f(one)./sum(C(one,:),2);
        tgt      = [y; tgt'];
        tgt_label = 'Target mean';
        if strcmpi(options.scenario.uncertainty,'hard')
            tgt_label = 'Exact target';
        end
    end
%%gda
end

nfigs  = ceil(length(varnames)/( nplots(1)*nplots(2)) );
nplots = repmat(nplots,nfigs,1);
for j=1:size(nplots,1),
    nbofplots(j)=nplots(j,1)*nplots(j,2);
end

ntotplots = sum(nbofplots);
if ntotplots<length(varnames),
    nfigplus = ceil((length(pplotvar)-ntotplots)/nbofplots(end));
    lastrow=nplots(end,:);
    lastrow=repmat(lastrow,nfigplus,1);
    nplots = [nplots;lastrow];
    nbofplots = [nbofplots repmat(nbofplots(end),1,nfigplus)];
end

if ndraws > 1    
    frcsts_ = nan(length(time),nvar,ndraws);
    for kk = 1 : ndraws
        frcsts_(:,:,kk) = [y; frcsts(:,:,kk)];
    end
    frcsts = frcsts_;
    
    sort_idx   = round((0.5 + [-conf_sig, conf_sig, 0]/2) * ndraws);
    
    if trasf_yes ==1 
        frcsts_ = nan(size(frcsts));
        for var = 1 : nvar 
            if order_trasform(var) == 1 % period over period
                frcsts_(2:end,var,:) = (frcsts(2:end,var,:) - frcsts(1:end-1,var,:)) ;
                
            elseif order_trasform(var) == 100 % percentage period over period 
                frcsts_(2:end,var,:) = 100*(frcsts(2:end,var,:) - frcsts(1:end-1,var,:));
                
            elseif order_trasform(var) == 12 % percentage 12 periord over 12 period (year over year % change f  or monthly data)
                frcsts_(13:end,var,:) = 100*(frcsts(13:end,var,:) - frcsts(1:end-12,var,:));
     
            elseif order_trasform(var) == 4 % percentage 4 periord over 4 period (year over year % change for quarterly data)
                frcsts_(5:end,var,:) = 100*(frcsts(5:end,var,:) - frcsts(1:end-4,var,:));
            
            elseif order_trasform(var) == 400
                frcsts_(2:end,var,:) = 400*(frcsts(2:end,var,:) - frcsts(1:end-1,var,:));
            
            elseif order_trasform(var) == 1200
                frcsts_(2:end,var,:) = 1200*(frcsts(2:end,var,:) - frcsts(1:end-1,var,:));
            
            else
                frcsts_(:,var,:) = frcsts(:,var,:);
            end
        end
        frcsts = frcsts_;
    end    
    
    frcsts_sort   = sort(frcsts,3);
    if sort_idx(1) == 0
        sort_idx(1) = 1;
        warning('Bands not reliable. You have too few draws.')
    end
    irf_Median = squeeze(frcsts_sort(:, :, sort_idx(3) ));
    irf_low    = squeeze(frcsts_sort(:, :, sort_idx(1) ));
    irf_up     = squeeze(frcsts_sort(:, :, sort_idx(2) ));
    if  add_multiple_bands_yes == 1
        irf_low_low  = squeeze(frcsts_sort(:, :, sort_idx_2(1) ));
        irf_up_up    = squeeze(frcsts_sort(:, :, sort_idx_2(2) ));
    end
    
else
    irf_Median = [y; frcsts];

    if trasf_yes ==1 
        irf_Median_ = nan(size(irf_Median));
        for var = 1 : nvar 
            if order_trasform(var) == 1
                irf_Median_(2:end,var) = diff(irf_Median(:,var)) ;
                
            elseif order_trasform(var) == 12 % percentage 12 periord over 12 period (year over year % change f  or monthly data)
                irf_Median_(13:end,var) = 100*(irf_Median(13:end,var) - irf_Median(1:end-12,var));
                
            elseif order_trasform(var) == 4 % percentage 4 periord over 4 period (year over year % change for quarterly data)
                irf_Median_(5:end,var) = 100*(irf_Median(5:end,var) - irf_Median(1:end-4,var));
            
            elseif order_trasform(var) == 100                
                irf_Median_(2:end,var) = 100*diff(irf_Median(:,var));
                
            elseif order_trasform(var) == 400
                irf_Median_(2:end,var) = 400*diff(irf_Median(:,var));
                
            elseif order_trasform(var) == 1200
                irf_Median_(2:end,var) = 1200*diff(irf_Median(:,var));
                
            else
                irf_Median_(:,var) = irf_Median(:,var);
            end
        end
        irf_Median = irf_Median_;
    end
    irf_low    = irf_Median;
    irf_up     = irf_Median;
    
    
end

if trasf_yes ==1 && add_frcst_yes ==1
    add_frcst_ = nan(size(add_frcst));
    for var = 1 : nvar
        if order_trasform(var) == 1
            add_frcst_(2:end,var) = diff(add_frcst(:,var)) ;
            
        elseif order_trasform(var) == 12 % percentage 12 periord over 12 period (year over year % change for monthly data)
            add_frcst_(13:end,var) = 100*(add_frcst(13:end,var) - add_frcst(1:end-12,var));
            
%         elseif order_trasform(var) == 4 % percentage 4 periord over 4 period (year over year % change for quarterly data)
%             add_frcst_(4:end,var) = 100*(add_frcst(4:end,var) - add_frcst(1:end-3,var));

        elseif order_trasform(var) == 4 % percentage 4 periord over 4 period (m over m % change for monthly data or YoY for Q data)
            add_frcst_(5:end,var) = 100*(add_frcst(5:end,var) - add_frcst(1:end-4,var));
            
        elseif order_trasform(var) == 100
            add_frcst_(2:end,var) = 100*diff(add_frcst(:,var));
            
        elseif order_trasform(var) == 400
            add_frcst_(2:end,var) = 400*diff(add_frcst(:,var));
            
        elseif order_trasform(var) == 1200
            add_frcst_(2:end,var) = 1200*diff(add_frcst(:,var));
            
        else
            add_frcst_(:,var) = add_frcst(:,var);
        end
    end
    add_frcst = add_frcst_;
end


%%gda
if scenario_yes == 1
    % the targets go through the same order_transform codes
    if trasf_yes == 1
        codes = [1 100 12 4 400 1200];  lags = [1 1 12 4 1 1];  mults = [1 100 100 100 400 1200];
        for var = 1 : nvar
            k = find(codes == order_trasform(var));
            if ~isempty(k)
                tgt(:,var) = [nan(lags(k),1); mults(k)*(tgt(lags(k)+1:end,var) - tgt(1:end-lags(k),var))];
            end
        end
    end
    tgt(1:length(y),:) = NaN;
    join = isfinite(tgt(length(y)+1,:));    % a target at period 1 is joined to the last observation
    tgt(length(y),join) = irf_Median(length(y),join);
end
%%gda
jplot = 0;
jfig  = 0;

index_time = time_start:time_end;

for var= 1: nvar
    
    if jplot==0,
%%gda
        figure('name',fig_name);    % was: figure('name',['Forecasts'] );
%%gda
        jfig=jfig+1;
    end
    
    jplot=jplot+1;
    subplot(nplots(1),nplots(2),jplot)    
    
    
%%gda
    if scenario_yes == 1
        % forecast rows start at the last observation
        jf = max(time_start,length(y)) : length(time);
        xf = time(jf);
        ih = time_start : length(y);
        shade_ = @(lo,up,col) patch('XData',[xf; flipud(xf)],'YData',[lo; flipud(up)],'FaceColor',col, ...
            'FaceAlpha',.9,'EdgeColor','none','HandleVisibility','off');
        hold on
        if add_multiple_bands_yes == 1
            shade_(irf_low_low(jf,var),irf_up_up(jf,var),[.80 .85 .97]);
            shade_(irf_low(jf,var),irf_up(jf,var),[.50 .60 .90]);
        else
            shade_(irf_low(jf,var),irf_up(jf,var),[.80 .85 .97]);
        end
        plot(time(ih),irf_Median(ih,var),'k','LineWidth',1.5,'DisplayName','History');
        if add_frcst_yes == 1
            plot(xf,add_frcst(jf,var),':','Color',[.45 .45 .45],'LineWidth',1.5,'DisplayName',add_label);
        end
        plot(xf,irf_Median(jf,var),'b-','LineWidth',2,'DisplayName','Median forecast');
        t   = tgt(jf,var);
        fin = isfinite(t);
        iso = fin & ~[false; fin(1:end-1)] & ~[fin(2:end); false];    % isolated targets are dots
        if any(fin & ~iso)
            plot(xf,t,'r-','LineWidth',1.2,'DisplayName',[tgt_label ' (each month)']);
        end
        if any(iso)
            plot(xf(iso),t(iso),'ro','LineStyle','none','MarkerFaceColor','r','MarkerSize',5, ...
                'DisplayName',[tgt_label ' (marked dates)']);
        end
        title(varnames{var})
        axis tight
        if jplot == 1 || any(fin)
            legend('show','Location','best','FontSize',8);
        end
        set(gcf,'position' ,[80 80 1400 760])
    else
%%gda
    if add_multiple_bands_yes == 1
        
        h = area([time(index_time)],[irf_low_low(index_time,var),...
            irf_low(index_time,var) - irf_low_low(index_time,var),...
            irf_up(index_time,var) - irf_low(index_time,var),...
            irf_up_up(index_time,var) - irf_up(index_time,var)]);%,'FaceColor',[.85 .85 .85]);
        set(h(4),'FaceColor',[.95 .95 .95])
        set(h(3),'FaceColor',[.85 .85 .85])
        set(h(2),'FaceColor',[.95 .95 .95])
        set(h(1),'FaceColor',[1 1 1])
        set(h,'linestyle','none')
        hold on
        min_=min(irf_low_low(index_time,var));
        max_=max(irf_up_up(index_time,var));
    else
        h = area([time(index_time)],[irf_low(index_time,var),...
            irf_up(index_time,var) - irf_low(index_time,var)]);%,'FaceColor',[.85 .85 .85]);
        set(h(2),'FaceColor',[.85 .85 .85])
        set(h(1),'FaceColor',[1 1 1])
        set(h,'linestyle','none')
        min_=min(irf_low(index_time,var));
        max_=max(irf_up(index_time,var));
    end
    hold on;
    plot(time(index_time),irf_Median(index_time,var),'k');
    hold on;
    plot(time(index_time),zeros(length(index_time),1),'k');
    hold on;
    axis tight
    ylim([min_,max_]);
    if add_frcst_yes == 1
        plot(time(index_time),add_frcst(index_time,var),'b','LineWidth',2);
        ylim([min([min_,min(add_frcst(index_time,var))]) , max([max_,max(add_frcst(index_time,var))])  ] );
    end
    title(varnames{var})
    set(gcf,'position' ,[50 50 800 650])
%%gda
    end
%%gda
    if jplot==nbofplots(jfig) || var==length(varnames)
        %legend(legenda);
        if savefig_yes == 1 
            STR_RECAP = [ fnam_dir '/' fnam_suffix  '_' int2str(jfig)];
            saveas(gcf,STR_RECAP,'fig');
            if strcmp(version('-release'),'2022b') == 0
                saveas(gcf,STR_RECAP,'eps');
                % savefigure_pdf(STR_RECAP);
                savefigure_pdf([STR_RECAP '.pdf']);
            end
        end
        jplot=0;
    end
end
