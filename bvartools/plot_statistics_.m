function plot_statistics_(y,options)

[T,ny] = size(y);
TT     = [1:1:T];

dirname = './';
tag = '';
savefig_yes = 0;
set_x_ticks = 0;
step_plot   = 20;
typek        = 1;
types        = 1;
n1 = 2 ; n2 = 2;
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
        savefig_yes = 1;
        dirname = options.dirname;
    end
    if isfield(options,'tag')
        savefig_yes = 1;
        tag = options.tag;
    end
    if isfield(options,'typek')
        typek = options.typek;
    end
    if isfield(options,'types')
        types = options.types;
    end
    if isfield(options,'TT')
        TT = options.TT;
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
    if isfield(options,'nplots') ==1;
        nplots = options.nplots;
        n1 = nplots(1);
        n2 = nplots(2);
    end
end


for ss = 1 : ny
    figure('Name',varnames{ss})
    subplot(n1,n2,1)
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
    subplot(n1,n2,2)
    probplot(y(:,ss))
    hold on
    subplot(n1,n2,3)
    h2 = histogram(y (:,ss),Kappa);
    
    [krt,krt_normal]  = kurtosis_(y(:,ss));
    [skw,skw_normal]  = skewness_(y(:,ss));
    ki = krt(typek) - krt_normal(typek);
    si = skw(types) - skw_normal(types);
    h2.Normalization = 'probability';%     h2.BinWidth = 0.25;
    title('Histogram')
    axis tight;
    xlabel(['ex k = ' num2str(ki,'%.1f') '; sk = ' num2str(si,'%.1f') ';'])
    if savefig_yes ==1
        STR_RECAP = [dirname '/stats_' tag varnames{ss}];
        savefigure_pdf([STR_RECAP '.pdf']);
        saveas(gcf,STR_RECAP,'fig');
    end
end
