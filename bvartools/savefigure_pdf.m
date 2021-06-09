function [ figformat ] = savefigure_pdf( name,varargin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018

% SAVEFIGURE_PDF: this function saves figures in .pdf format
% name: Name of the figure, example: 'figure1'
% type: Format of the figure, example: 'fig','espc'
% If no input is provided, figure saved twice with fig and espc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,name,'-dpdf','-r0');
figformat=1;
end

