function [same, maxabs, info] = bitcompare(a, b)
%BITCOMPARE  Compare two values for exact equality, recursing into cells.
%
%   The extracted functions are supposed to be the reference's own code, so
%   the bar here is equality, not tolerance. Any difference at all means the
%   hoist changed something. MAXABS is reported anyway, because "differs by
%   3e-16" and "differs by 0.4" call for very different investigations.
%
%   NaN is treated as equal to NaN: the reference produces NaN legitimately
%   where a covariance ratio is degenerate, and a gate that called those
%   unequal would fail on correct output.

if iscell(a) && iscell(b)
    if ~isequal(size(a), size(b))
        same = false; maxabs = Inf;
        info = sprintf('cell size %s vs %s', mat2str(size(a)), mat2str(size(b)));
        return
    end
    same = true; maxabs = 0;
    for i = 1:numel(a)
        [element_same, element_max] = bitcompare(a{i}, b{i});
        same = same && element_same;
        maxabs = max(maxabs, element_max);
    end
    info = sprintf('cell[%d]', numel(a));

elseif isnumeric(a) && isnumeric(b)
    if ~isequal(size(a), size(b))
        same = false; maxabs = Inf;
        info = sprintf('size %s vs %s', mat2str(size(a)), mat2str(size(b)));
        return
    end
    same = isequaln(a, b);
    difference = abs(a - b);
    difference(isnan(difference)) = 0;
    if isempty(difference)
        maxabs = 0;
        info = 'empty';
    else
        maxabs = max(difference(:));
        info = sprintf('numeric%s', mat2str(size(a)));
    end

else
    same = isequaln(a, b);
    maxabs = NaN;
    info = sprintf('%s vs %s', class(a), class(b));
end

end
