function [share, finalCond] = hypothesisHump(x, maxIdx, dumStrong)
[~, idx] = max(x, [], 2);
weakCond = idx > maxIdx & idx < size(x,2);

finalCond = weakCond;
if dumStrong == 1
    for dd = 1:numel(idx)
        left = diff(x(dd,1:idx(dd)));
        right = diff(x(dd,idx(dd):end));
        aux = all(left > 0) & all(right < 0);
        finalCond(dd) = aux & finalCond(dd);
    end
end

share = sum(finalCond) / numel(idx) * 100;

end