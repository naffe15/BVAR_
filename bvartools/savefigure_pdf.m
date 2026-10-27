function [ figformat ] = savefigure_pdf( name,varargin )
%SAVEFIGURE_PDF Saves the current figure as a PDF cropped to the figure size.
%
%   SAVEFIGURE_PDF(name)
%
%   Inputs
%     name   file name, with or without path and extension
%
%   Output
%     figformat   1
%
%   Further inputs are accepted and ignored.
%
%   See also PLOT_IRFS_, PLOT_FRCST_.

%   Filippo Ferroni, 6/1/2015.

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,name,'-dpdf','-r0');
figformat=1;
end

