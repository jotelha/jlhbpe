% functionality of char array functions below is to be replaced by
% functions working on cell arrays, such as:
function out = cprintf(format,c) 
    out = cellfun(@(s) sprintf(format,s),c,'UniformOutput',false);
end