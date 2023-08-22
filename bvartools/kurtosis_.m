function [krt,krt_normal]  = kurtosis_(y,flags_)

[~,ny] = size(y);
if nargin<2
    flags_ =ones(ny,4);
end

% This function computes the kurtosis using three different
% estimators which are robust to outliers
% =====================
% 1. k1 = ratio of fourth (centered) moment / square of the second
% (centered) moment
% NORMAL: k1 = 3
% =====================
% 2. k2 = ((p(0.875) - p(0.625)) + (p(0.375)-p(0.125)))/(p(0.75)-p(0.25))
% where p(j) is the j-percentile
% See Moors(1988) A quantile alternative for kurtosis, The Statistician, 37
% (1988), pp. 25-32
% NORMAL: k2 = 1.23;
% =====================
% 3. k3 = (U(0.05)-L(0.05)) / (U(0.5)-L(0.5))
% See Hoog (1976) More light on the kurtosis and related statistics Journal of the American
% Statistical Association, 67 (1972), pp. 422-424
% NORMAL: k3 = 2.95;
% =====================
% 4. k4 = (F^-1(0.975)-F^-1(0.025)) / (F^-1(0.75)-F^-1(0.25))
% NORMAL: k4= 2.91;
% See Crow and Siddiqui, 1967, Robust estimation of location, Journal of the
% American Statistical Association, 62 (1967), pp. 353-389

krt = nan(ny,4);
%p1 = [0.15    0.2500    0.3750    0.5000    0.6250    0.7500    0.85];
p1 = [12.5 : 12.5 : 87.5]/100;
p2 = [0.025 0.25 0.75 0.975];

krt_normal(1) = 3;
krt_normal(2) = 1.23; %((norminv(p1(7))-norminv(p1(5))) + (norminv(p1(3))-norminv(p1(1))) )/ (norminv(p1(6))-norminv(p1(2)));
krt_normal(3) = 2.59;
krt_normal(4) = 2.90; %(norminv(0.975)-norminv(0.025)) / (norminv(0.75)-norminv(0.25));


for j = 1 : ny
    y_ = y(:,j);
    % standard
    if (flags_ (j,1)) == 1
        krt(j,1) = kurtosis(y(:,j),0);
    end    
    if any(flags_ (j,2:end)) == 1
        y_(isnan(y_))= [];
        T = size(y_,1);
        ysort = sort(y_);
    end
    % Moors
    if (flags_ (j,2)) == 1
        emp_(1) = ysort(floor(p1(1)*T));
        emp_(2) = ysort(ceil(p1(2)*T));
        emp_(3) = ysort(ceil(p1(3)*T));
        emp_(4) = ysort(ceil(p1(4)*T));
        emp_(5) = ysort(floor(p1(5)*T));
        emp_(6) = ysort(floor(p1(6)*T));
        emp_(7) = ysort(ceil(p1(7)*T));
        krt(j,2) = ((emp_(7)-emp_(5)) + (emp_(3)-emp_(1)) )/ (emp_(6)-emp_(2));
    end
    % Hogg
    if (flags_ (j,3)) == 1
        % krt(j,3) = hogg(y(:,j));
        u1 = floor(T*0.95);
        l1 = ceil(T*0.05);
        u2 = round(T*0.5);
        l2 = round(T*0.5);
        U1 = mean(ysort(u1:end));
        L1 = mean(ysort(1:l1));
        U2 = mean(ysort(u2:end));
        L2 = mean(ysort(1:l2));
        krt(j,3) = (U1-L1) / (U2-L2);
    end
    % Crow and Siddiqui
    if (flags_ (j,4)) == 1
        emp2_(1) = ysort(ceil(p2(1)*T));
        emp2_(2) = ysort(round(p2(2)*T));
        emp2_(3) = ysort(round(p2(3)*T));
        emp2_(4) = ysort(floor(p2(4)*T));
        krt(j,4) = (emp2_(4)- emp2_(1)) / (emp2_(3) - emp2_(2));
    end
end



%     function [K] = hogg(y)
%         y(isnan(y))= [];
%         T = size(y,1);
%         u1 = floor(T*0.95);
%         l1 = ceil(T*0.05);
%         u2 = round(T*0.5);
%         l2 = round(T*0.5);
%         ysort = sort(y);
%         U1 = mean(ysort(u1:end));
%         L1 = mean(ysort(1:l1));
%         U2 = mean(ysort(u2:end));
%         L2 = mean(ysort(1:l2));
%         K = (U1-L1) / (U2-L2);
%     end
%
%     function [K] = moors(y)
%         y(isnan(y))= [];
%         p1 = [12.5 : 12.5 : 87.5]/100;
%         emp_   = prctile(y,p1);
%         K = ((emp_(7)-emp_(5)) + (emp_(3)-emp_(1)) )/ (emp_(6)-emp_(2));
%     end

end