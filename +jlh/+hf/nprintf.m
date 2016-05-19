function out = nprintf(format,n) 
    out = arrayfun(@(i) sprintf(format,i),(1:n)','UniformOutput',false);
end
