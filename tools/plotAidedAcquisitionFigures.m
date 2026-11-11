function [figures,files] = plotAidedAcquisitionFigures(maps,rec,cfg,name)
% 'separate': two 1x3 figures; 'six': one 2x3 figure with shared column captions.
% Panels are B1I, B1C, joint. Data/peak values are not modified.
required={'B1ICodeTimeMs','B1CCodeTimeMs','jointCodeTimeMs'};
assert(all(isfield(rec,required)), ...
    'This cache predates exact code-time peak coordinates; rerun acquisition.');
if cfg.saveFigures && ~isfolder(cfg.resultsDir), mkdir(cfg.resultsDir); end
dpi=300; if isfield(cfg,'figureResolution'), dpi=cfg.figureResolution; end
assert(isscalar(dpi) && isfinite(dpi) && dpi>=100,'Invalid figure resolution.');
figureLayout='separate';
if isfield(cfg,'figureLayout'), figureLayout=lower(string(cfg.figureLayout)); end
assert(isscalar(string(figureLayout)) && ismember(string(figureLayout),["separate","six"]), ...
    'figureLayout must be separate or six.');
sixPanels=strcmp(figureLayout,'six');
phases={maps.phaseAxisI,maps.phaseAxisC,maps.phaseAxisI};
dopplers={maps.fdAxisI,maps.fdAxisC,maps.fdAxisC};
surfaces={maps.surfaceI,maps.surfaceC,maps.surfaceJ};
peaks=[rec.B1IDopplerHz rec.B1ICodeTimeMs rec.B1IRatio; ...
    rec.B1CDopplerHz rec.B1CCodeTimeMs rec.B1CRatio; ...
    rec.jointB1CDopplerHz rec.jointCodeTimeMs rec.jointScore];
highest=max(peaks(:,3));
for k=1:3, highest=max(highest,max(double(surfaces{k}(:)))); end
assert(isfinite(highest) && highest>0,'Invalid metric surface.');
step=10^floor(log10(highest))/2;
zTop=ceil(highest*1.08/step)*step;
if sixPanels, tags={'Six_combined'}; else, tags={'3D_combined','Side_combined'}; end
cameras=[-37.5 30; -45 10];
captions={'(a) B1I','(b) B1C','(c) Joint acquisition'};
figures=gobjects(1,numel(tags));
files=repmat(struct('branch','','png','','fig',''),1,numel(tags));
for figureIndex=1:numel(tags)
    height=11.5; if sixPanels, height=22; end
    fig=figure('Name',[name '_' tags{figureIndex}],'Color','w','Units','centimeters', ...
        'Position',[1+.5*(figureIndex-1) 2+.5*(figureIndex-1) 28 height], ...
        'Visible',cfg.figureVisible,'InvertHardcopy','off','NumberTitle','off');
    figures(figureIndex)=fig;
    views=figureIndex; if sixPanels, views=1:2; end
    for viewIndex=views
    labelY=.115; labelHeight=.040;
    if sixPanels
        labelHeight=.022;
        if viewIndex==1
            panelY=.615; labelY=.548;
        else
            panelY=.125; labelY=.052;
        end
    else
        layout=tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');
        layout.OuterPosition=[.005 .19 .99 .81];
    end
    for k=1:3
        if sixPanels
            ax=axes(fig,'Units','normalized','Position',[(k-1)/3+.045 panelY .255 .345]);
        else
            ax=nexttile(layout);
        end
        phase=double(phases{k}(:)); fd=double(dopplers{k}(:));
        z=double(surfaces{k});
        % The reference uses Doppler on X and code delay on Y.
        if isscalar(fd)
            plot3(ax,repmat(fd,size(phase)),phase,z(:),'Color',[0 0 .7],'LineWidth',.6);
        else
            surf(ax,fd,phase,z.','FaceColor','interp', ...
                'EdgeColor',[0 .015 .08],'EdgeAlpha',.10,'LineWidth',.12, ...
                'FaceLighting','none');
        end
        view(ax,cameras(viewIndex,:)); grid(ax,'on');
        ax.FontName='Times New Roman'; ax.FontSize=9;
        ax.LineWidth=.65; ax.Color='w'; ax.Box='off';
        ax.GridColor=[0 0 0]; ax.GridLineStyle=':'; ax.GridAlpha=.55;
        ax.XColor='k'; ax.YColor='k'; ax.ZColor='k';
        ax.TickDir='in'; ax.Projection='orthographic';
        zlabel(ax,'Normalized metric','FontName','Times New Roman','FontSize',11);
        % Common metric/color scales across all methods and both views.
        phaseLimits=[0 1]; if k==2, phaseLimits=[0 10]; end
        ylim(ax,phaseLimits); zlim(ax,[0 zTop]); clim(ax,[0 zTop]);
        requested=double(cfg.dopplerHz(:));
        if k==1, requested=requested*(1561.098/1575.42); end
        fdLimits=[min([requested;fd;peaks(k,1)]) max([requested;fd;peaks(k,1)])];
        if diff(fdLimits)==0, fdLimits=fdLimits+[-.5 .5]; end
        xlim(ax,fdLimits);
        if viewIndex==2
            ax.YDir='reverse';
            % Equal base edges and a 45-degree oblique orthographic camera
            % give equal screen-projected Doppler/code-axis lengths.
            pbaspect(ax,[1 1 1.25]);
            view(ax,cameras(viewIndex,:));
        end
        ax.XTick=linspace(fdLimits(1),fdLimits(2),5);
        ax.YTick=linspace(phaseLimits(1),phaseLimits(2),6);
        xtickformat(ax,'%.0f'); ytickformat(ax,'%.1f');
        ax.XTickLabelMode='auto';
        if viewIndex==2
            ax.XTickLabelRotation=0; ax.YTickLabelRotation=0;
            % At the front corner the two endpoint labels meet. Move the
            % Doppler labels outward rather than hiding any of them.
            if isprop(ax.XAxis,'TickLabelGapOffset'), ax.XAxis.TickLabelGapOffset=8; end
            if isprop(ax.YAxis,'TickLabelGapOffset'), ax.YAxis.TickLabelGapOffset=0; end
        end
        ax.XAxis.Exponent=0; ax.YAxis.Exponent=0; ax.ZAxis.Exponent=0;
        colormap(ax,jet(256));
        ax.Tag=sprintf('AcquisitionPanel_%d_%d',viewIndex,k);
        marker=markAcquisitionPeak(ax,peaks(k,:),false);
        % X/Y order differs from the original single-panel helper.
        rows=marker.DataTipTemplate.DataTipRows;
        rows(1).Label='Doppler (Hz)'; rows(1).Format='%.3f';
        rows(2).Label='Code delay (ms)'; rows(2).Format='%.6f';
        marker.DataTipTemplate.DataTipRows=rows;
        marker.MarkerSize=4.5; marker.LineWidth=.7;
        % Fixed 2-D label positions avoid overlapping 3-D labels between
        % neighbouring panels, especially in the near-horizontal view.
        annotation(fig,'textbox',[(k-1)/3+.012 labelY .17 labelHeight], ...
            'String','Code delay (ms)','FontName','Times New Roman','FontSize',11, ...
            'HorizontalAlignment','center','EdgeColor','none','Margin',0);
        annotation(fig,'textbox',[(k-1)/3+.19 labelY .135 labelHeight], ...
            'String','Doppler (Hz)','FontName','Times New Roman','FontSize',11, ...
            'HorizontalAlignment','center','EdgeColor','none','Margin',0);
    end
    end
    % One caption per method column, shared by the upper/lower panels.
    captionY=.018; captionHeight=.045;
    if sixPanels, captionY=.008; captionHeight=.026; end
    for k=1:3
        annotation(fig,'textbox',[(k-1)/3+.015 captionY 1/3-.03 captionHeight], ...
            'String',captions{k},'FontName','Times New Roman','FontSize',12, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'EdgeColor','none','Margin',0,'Interpreter','none','Tag','AcquisitionColumnCaption');
    end
    drawnow;
    files(figureIndex).branch=tags{figureIndex};
    if cfg.saveFigures
        files(figureIndex).png=fullfile(cfg.resultsDir,[name '_' tags{figureIndex} '.png']);
        files(figureIndex).fig=fullfile(cfg.resultsDir,[name '_' tags{figureIndex} '.fig']);
        exportgraphics(fig,files(figureIndex).png,'Resolution',dpi,'BackgroundColor','white');
        savefig(fig,files(figureIndex).fig);
    end
end
end
