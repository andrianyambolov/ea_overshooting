function [share, cond] = hypothesisMax(x, maxIdx)
[~, idx] = max(x, [], 2);

cond = (idx <= maxIdx);
share = sum(cond) / numel(idx) * 100;
end

