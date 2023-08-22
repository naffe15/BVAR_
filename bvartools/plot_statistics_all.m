function plot_statistics_all(y,options)

[T,ny] = size(y);
TT     = [1:1:T];

dirname = './';
tag = '';
set_x_ticks = 0;
step_plot   = 20;
typek = 4;
types = 4;
save_figure = 0;

% declaring the names for the observable variables
for v = 1 : ny
    eval(['varnames{'   num2str(v) '} =  ''Var' num2str(v) ''';'])
end
Kappa = ceil (1+ log2(length(y)) );

if nargin > 1
    if isfield(options,'Kappa')
        Kappa = options.Kappa;
    end
    if isfield(options,'dirname')
        dirname = options.dirname;
        save_figure = 1;
    end
    if isfield(options,'saveas_dir')
        dirname = options.saveas_dir;
        save_figure = 1;
    end
    if isfield(options,'tag')
        tag = options.tag;
        save_figure = 1;
    end
    if isfield(options,'TT')
        TT = options.TT;
    end
    if isfield(options,'typek')
        typek = options.typek;
    end
    if isfield(options,'types')
        types = options.types;
    end
    if isfield(options,'set_x_ticks') && options.set_x_ticks ==1
        set_x_ticks = 1;
        time = datestr(TT);
    end
    
    if isfield(options,'step_plot')
        step_plot = options.step_plot;
    end
    if isfield(options,'varnames')==1
        varnames = options.varnames;
        if length(varnames) ~= ny
            error('wrong number of varnames')
        end
    end
end

figure,
for ss = 1 : ny
    
    %figure('Name',varnames{ss})
    subplot(2,ny,ss)
    plot(TT,y(:,ss))
    axis tight
    grid on
    hold on
    title(varnames{ss})
    if set_x_ticks == 1
        set(gca,'Xtick',TT(1:step_plot:end))
        tmp_str= time(1:step_plot:end,:);
        set(gca,'Xticklabel',tmp_str)
    end
    subplot(2,ny,ny + ss)
    probplot(y(:,ss))
    hold on
    axis tight
    title('')
    [krt,krt_normal]  = kurtosis_(y(:,ss));
    [skw,skw_normal]  = skewness_(y(:,ss));
    ki = krt(typek) - krt_normal(typek);
    si = skw(types) - skw_normal(types);
    xlabel(['ex k = ' num2str(ki,'%.1f') '; sk = ' num2str(si,'%.1f') ';'])
end
set(gcf,'position' ,[50 50 1100 550])
if save_figure
    STR_RECAP = [dirname '/stats_' tag];
    savefigure_pdf([STR_RECAP '.pdf']);
    saveas(gcf,STR_RECAP,'fig');
end