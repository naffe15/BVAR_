function [skw,skw_normal]  = skewness_(y,flags_)

[~,ny] = size(y);
if nargin<2
    flags_ =ones(ny,4);
end

% This function computes the skewness using three different
% estimators which are robust to outliers
% =====================
% 1. s1 = ratio of third (centered) moment / cube of the standard
% deviation
% NORMAL: k1 = 0;
% =====================
% 2. s2 = ((p(0.75) + p(0.25) - 2*p(0.50))/(p(0.75)-p(0.25))
% where p(j) is the j-percentile
% See A.L. Bowley Elements of Statistics, Scribner's, New York (1920)
% NORMAL: s2 = 0;
% =====================
% 3. s3 = (mean-p(0.5)) / E|y(t)-p(0.50)|
% See R.A. Groeneveld, G. Meeden, Measuring skewness and kurtosis The
% Statistician, 33 (1984), pp. 391-399 
% NORMAL: s3 = 0;
% =====================
% 4. s4 = (mean-p(0.5)) / sigma
% NORMAL: s4 = 0;
% See M.G. Kendall, A. Stuart, The Advanced Theory of Statistics, vol. 1,
% Griffin, London (1977) 

skw = nan(ny,4);
p1 = [0.25 0.5 0.75];
skw_normal = zeros(4,1);

for j = 1 : ny
    y_ = y(:,j);
        
    % standard
    if flags_ (j,1) == 1
        skw(j,1) = skewness(y(:,j),0);
    end
    if any(flags_ (j,2:end)) == 1
        y_(isnan(y_))= [];
        T = size(y_,1);
        ysort = sort(y_);
        emp_(1) = ysort(ceil(p1(1)*T));
        emp_(2) = ysort(ceil(p1(2)*T));
        emp_(3) = ysort(ceil(p1(3)*T));
    end
    % Bowley
    if (flags_ (j,2)) == 1
        skw(j,2) = (emp_(3)  + emp_(1) - 2*emp_(2)) / (emp_(3)-emp_(1));
    end
    % Groeneveld-Meeden
    if (flags_ (j,3)) == 1
        meany  = mean(y_);
        absdev   = 1/length(y_) * sum(abs(y_-emp_(2)));
        skw(j,3) = (meany-emp_(2)) / absdev;
    end    
    % Kendall-Stuart
    if (flags_ (j,4)) == 1
        meany  = mean(y_);
        stdevy   = std(y_);
        skw(j,4) = (meany-emp_(2)) / stdevy;
    end
end

end