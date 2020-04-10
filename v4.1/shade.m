function []=shade(start,finish,colorstr,up,low);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filippo Ferroni, 6/1/2015
% Revised, 2/15/2017
% Revised, 3/21/2018

% function []=shade(start,finish,colorstr);
%
%  start and finish are Nx1 vectors of starting and ending years.
%  The function shades between the start and finish pairs using colorstr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('colorstr'); colorstr='y'; end;  % default is yellow
curax=axis;
y=[curax(3)+low up*curax(4) up*curax(4) curax(3)+low];
hold on;
for i=1:length(start);
  x=[start(i) start(i) finish(i) finish(i)];
  fill(x,y,colorstr);
end;
  
% Now, prevent the shading from covering up the lines in the plot.  
h = findobj(gca,'Type','line');
set(h,'EraseMode','xor');

h = findobj(gca,'Type','patch');
set(h,'EdgeColor','none');

% This last one makes the tick marks visible
set(gca, 'Layer', 'top')
