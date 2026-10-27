function []=shade(start,finish,colorstr,up,low);
%SHADE Shades periods (e.g. recessions) on the current plot.
%
%   SHADE(start,finish,colorstr,up,low)
%
%   Inputs
%     start, finish   N x 1   first and last date of each shaded period
%     colorstr        colour (default 'y')
%     up              top of the band, as a fraction of the upper y limit
%     low             offset added to the lower y limit
%
%   See also PLOT_FRCST_, VLINE.

%   Filippo Ferroni, 6/1/2015.

if ~exist('colorstr'); colorstr='y'; end;  % default is yellow
curax=axis;
y=[curax(3)+low up*curax(4) up*curax(4) curax(3)+low];
hold on;
for i=1:length(start);
  x=[start(i) start(i) finish(i) finish(i)];
  fill(x,y,colorstr);
end;
  
% Now, prevent the shading from covering up the lines in the plot.  
% h = findobj(gca,'Type','line');   % old: EraseMode was removed in MATLAB R2014b
% set(h,'EraseMode','xor');          % and errors there
ch  = get(gca,'Children');           % move the shaded areas behind everything else,
isp = strcmp(get(ch,'Type'),'patch');%   so lines plotted before shade stay visible
set(gca,'Children',[ch(~isp); ch(isp)]);

h = findobj(gca,'Type','patch');
set(h,'EdgeColor','none');

% This last one makes the tick marks visible
set(gca, 'Layer', 'top')
