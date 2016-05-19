function out = aprintf(format,a) 
    out = arrayfun(@(a) sprintf(format,a),a,'UniformOutput',false);
end
