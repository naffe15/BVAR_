function xq = m2q(x,options)
%M2Q Monthly to quarterly aggregation of one series.
%
%   xq = M2Q(x)
%   xq = M2Q(x,options)
%
%   Inputs
%     x        Tm x 1   monthly series, starting in the first month of a
%                       quarter
%     options  struct (optional); the presence of a field selects the
%              aggregation (default: mean of the three months, NaN ignored)
%       .median   median of the three months
%       .end      last month of the quarter
%       .start    first month of the quarter
%
%   Output
%     xq       floor(Tm/3) x 1   quarterly series

if size(x,1) < size(x,2)
    x = x';
end
if size(x,2) > 1
    error('x should be a column vector')
end


% start = 1;
Tm    = length(x);
Tq    = floor(Tm/3);
aggr  = 1; % mean


if nargin > 1
    if isfield(options,'median') ==1
        % take the median over the quarter
        aggr = 2;
    end
    if isfield(options,'end') ==1
        % take the end of the quarter
        aggr = 3;
    end
    if isfield(options,'start') ==1
        % take the start of the quarter
        aggr = 4;
    end
end

xq = nan(Tq,1);

switch aggr
    case 1 % mean
        for tt = 1 : Tq
            xq(tt,1) = nanmean(x(3*(tt-1)+1 : 3*tt,1));
        end
    case 2 % median
        for tt = 1 : Tq
            xq(tt,1) = nanmedian(x(3*(tt-1)+1 : 3*tt,1));
        end
    case 3 % end of q
        for tt = 1 : Tq
            xq(tt,1) = x(3*tt,1);
        end
    case 4 % start of q
        for tt = 1 : Tq            
            xq(tt,1) = x(3*(tt-1)+1 ,1);
        end
end


