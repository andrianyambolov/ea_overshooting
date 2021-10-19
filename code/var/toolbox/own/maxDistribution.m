function [share, cond] = maxDistribution(x)
[~, idx] = max(x, [], 2);

cond = [];
for jj = 1:size(x,2)
    aux = sum(idx == jj);
    cond = [cond; aux];
end
share = cumsum(cond) / numel(idx) * 100;
end

