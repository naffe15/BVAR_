function [x,Scale]=standard(y)
T=size(y,1);
my=repmat(nanmean(y),T,1);
sy=repmat(nanstd(y),T,1);
x=(y-my)./sy;
x(find(isnan(x)))=randn;
Scale = nanstd(y)';
%x=(y-kron(mean(y),ones(rows(y),1)))./kron(std(y),ones(rows(y),1));
