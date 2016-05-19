% uses a cell array of strings to make a diagonal cell matrix, like diag
% for numericak arrays
function out = strdiag(c) 
    import fp.iif
    out = arrayfun( ... 
        @(p) iif(  p == 0, @() '0', true, @() c{p} ), ...
        diag(1:numel(c)), 'UniformOutput',false);
end