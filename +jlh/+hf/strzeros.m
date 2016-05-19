function out = strzeros(n,m) 
    out = arrayfun( @(p) '0', zeros(n,m), 'UniformOutput',false);
end