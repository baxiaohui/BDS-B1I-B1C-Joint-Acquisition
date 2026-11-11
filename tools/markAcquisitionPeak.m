function [marker,tip] = markAcquisitionPeak(ax,peakXYZ,showTip)
%MARKACQUISITIONPEAK Mark the global fine-grid maximum on an acquisition plot.
% peakXYZ = [code-time delay (ms), branch Doppler (Hz), normalized score].
% The caller supplies the exact full-search winner; do not infer its
% coordinates from a downsampled / max-pooled display surface.
assert(isgraphics(ax,'axes'),'A valid acquisition axes is required.');
assert(isnumeric(peakXYZ) && isreal(peakXYZ) && numel(peakXYZ)==3 && ...
    all(isfinite(peakXYZ)),'Peak must be finite [delay_ms Doppler_Hz score].');
if nargin<3, showTip=true; end
% Idempotent: refreshing a figure must not accumulate markers or labels.
delete(findall(ax,'Tag','AcquisitionMaximumMarker'));
wasHeld=ishold(ax); hold(ax,'on');
marker=plot3(ax,peakXYZ(1),peakXYZ(2),peakXYZ(3),'o', ...
    'MarkerSize',7,'MarkerFaceColor',[.15 .15 .15],'MarkerEdgeColor','w', ...
    'LineWidth',1.2,'Tag','AcquisitionMaximumMarker','HandleVisibility','off');
if ~wasHeld, hold(ax,'off'); end
rows=marker.DataTipTemplate.DataTipRows;
rows(1).Label='Code delay (ms)'; rows(1).Format='%.6f';
rows(2).Label='Doppler (Hz)'; rows(2).Format='%.3f';
rows(3).Label='Peak'; rows(3).Format='%.3f';
marker.DataTipTemplate.DataTipRows=rows;
% Only a small margin is needed when the persistent popup is disabled.
limits=zlim(ax); span=max([peakXYZ(3)-limits(1),abs(peakXYZ(3))*.1,eps]);
margin=.05; if showTip, margin=.25; end
if limits(2)<peakXYZ(3)+margin*span
    zlim(ax,[limits(1),peakXYZ(3)+margin*span]);
end
tip=[];
if showTip, tip=datatip(marker,'DataIndex',1,'FontSize',9,'Interpreter','none'); end
end
