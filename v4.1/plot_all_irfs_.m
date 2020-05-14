function plot_all_irfs_(irfs,options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018

% input : irfs
% 1st dimension: variable
% 2nd dimension: horizon
% 3rd dimension: shock
% 4th dimension: draws
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvar    = size(irfs,1);
hor     = size(irfs,2);
nshocks = size(irfs,3);
ndraws  = size(irfs,4);
nplots  = [nvar nshocks];
savefig_yes = 0;
conf_sig    = 0.68;
normz       = 1;
add_irfs_yes = 0;
normz_yes   = 0;
add_multiple_bands_yes = 0;

if nargin <2
    disp('You did not provided names for shocks nor variables.')
    disp('I call them Var 1, Var 2, ... and Shck 1, ...')
    for v = 1 : nvar
        eval(['varnames{'   num2str(v) '} =  ''Var  ' num2str(v) ''';'])
    end
    for v = 1 : nshocks
        eval(['shocksnames{' num2str(v) '} =  ''Shck ' num2str(v) ''';'])
    end
else
    if isfield(options,'varnames') ==1
        varnames = options.varnames;
    end
    if isfield(options,'shocksnames') == 1
        shocksnames = options.shocksnames;
        if length(shocksnames) ~= nshocks
            error('There is a mismatch between the number of shocks and the names')
        end
    else
        for v = 1 : nshocks
            eval(['shocksnames{' num2str(v) '} =  ''Shck ' varnames{v} ''';'])
        end
    end
    if isfield(options,'normz') ==1 && options.normz==1,
        normz_yes = 1;
    end
    if isfield(options,'conf_sig') ==1;
        conf_sig = options.conf_sig;
    end
    if isfield(options,'nplots') ==1;
        nplots = options.nplots;
    end
    if isfield(options,'saveas_strng') ==1;
        savefig_yes = 1;
        % setting the names of the figure to save
        fnam_suffix = [ 'irfs_' options.saveas_strng ];
        fnam_dir    = '.';
    end
    if isfield(options,'saveas_dir') ==1;
        % setting the folder where to save the figure
        fnam_dir = options.saveas_dir;
        if exist(fnam_dir,'dir') == 0
            mkdir(fnam_dir)
        end
    end
    if isfield(options,'add_irfs') ==1
        % setting the folder where to save the figure
        add_irfs = options.add_irfs;
        add_irfs_yes = 1;
    end
    if isfield(options,'conf_sig_2') ==1
        add_multiple_bands_yes = 1;
        sort_idx_2   = round((0.5 + [-options.conf_sig_2, options.conf_sig_2, 0]/2) * ndraws);
    end

end



% nfigs  = ceil(length(varnames)/( nplots(1)*nplots(2)) );
% nplots = repmat(nplots,nfigs,1);
% for j=1:size(nplots,1),
%     nbofplots(j)=nplots(j,1)*nplots(j,2);
% end
% 
% ntotplots = sum(nbofplots);
% if ntotplots<length(varnames),
%     nfigplus = ceil((length(pplotvar)-ntotplots)/nbofplots(end));
%     lastrow=nplots(end,:);
%     lastrow=repmat(lastrow,nfigplus,1);
%     nplots = [nplots;lastrow];
%     nbofplots = [nbofplots repmat(nbofplots(end),1,nfigplus)];
% end
% 
% conf_sig = 0.68;
% sort_idx = round((0.5 + [-conf_sig, conf_sig, 0]/2) * options.K);
%
% sims_shock_down_conf = normz * sims_shock_sort(:, :, :,  sort_idx(1));
% sims_shock_up_conf   = normz * sims_shock_sort(:, :, :,  sort_idx(2));
% sims_shock_median    = normz * sims_shock_sort(:, :, :,  sort_idx(3));

if ndraws > 1
    sort_idx   = round((0.5 + [-conf_sig, conf_sig, 0]/2) * ndraws);
    irf_sort   = sort(irfs,4);
    if normz_yes == 1
        % normalize the IRF relative to a 100 bpt increase in the first
        % variable, first horizon, first shock
        normz = 1/ irf_sort(1, 1, 1,  sort_idx(3) );
    end
    if sort_idx(1) == 0,
        sort_idx(1) = 1;
        warning('IRF Bands not reliable. You have too few draws.')
    end
    irf_Median = normz * squeeze(irf_sort(:, :, :,  sort_idx(3) ));
    irf_low    = normz * squeeze(irf_sort(:, :, :,  sort_idx(1) ));
    irf_up     = normz * squeeze(irf_sort(:, :, :,  sort_idx(2) ));
    if  add_multiple_bands_yes == 1
        irf_low_low  = normz * squeeze(irf_sort(:, :, :,  sort_idx_2(1) ));
        irf_up_up    = normz * squeeze(irf_sort(:, :, :,  sort_idx_2(2) ));
    end

else
    irf_Median = normz * irfs;
    irf_low    = normz * irfs;
    irf_up     = normz * irfs;
end

jplot = 0;
% jfig  = 0;
figure('name',['All IRFs'] );
for sho = 1 : nshocks
    for var= 1: size(varnames,2)
        
        
        jplot=jplot+1;
        subplot(nplots(1),nplots(2),jplot)
        
        if add_multiple_bands_yes == 1            
            h = area([irf_low_low(var,:,sho)',...
                irf_low(var,:,sho)' - irf_low_low(var,:,sho)',...
                irf_up(var,:,sho)' - irf_low(var,:,sho)',...
                irf_up_up(var,:,sho)' - irf_up(var,:,sho)']);%,'FaceColor',[.85 .85 .85]);
            set(h(4),'FaceColor',[.95 .95 .95])
            set(h(3),'FaceColor',[.85 .85 .85])
            set(h(2),'FaceColor',[.95 .95 .95])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
            hold on
            
        else
            h = area([irf_low(var,:,sho)',...
                irf_up(var,:,sho)' - irf_low(var,:,sho)']);%,'FaceColor',[.85 .85 .85]);
            set(h(2),'FaceColor',[.85 .85 .85])
            set(h(1),'FaceColor',[1 1 1])
            set(h,'linestyle','none')
            hold on;
        end
        
        plot(irf_Median(var,:,sho),'k');
        if add_irfs_yes == 1
            plot(add_irfs(var,:,sho),'b','LineWidth',2);
        end
%         if add_multiple_bands_yes == 1
%             plot(irf_up_up(var,:,sho),'k:','LineWidth',1.2);
%             plot(irf_low_low(var,:,sho),'k:','LineWidth',1.2);
%         end
        hold on;
        plot(zeros(1,hor),'k')
        hold on;
        hold on
        axis tight
        if jplot <= nvar
            title(varnames{var})
        end
        if jplot == nvar*(sho-1) + 1
            ylabel(shocksnames{sho})
        end
        set(gcf,'position' ,[50 50 800 650])
    end
end
if savefig_yes == 1,
    STR_RECAP = [ fnam_dir '/' fnam_suffix ];
    saveas(gcf,STR_RECAP,'fig');
    saveas(gcf,STR_RECAP,'eps');
    savefigure_pdf([STR_RECAP '.pdf']);
end
