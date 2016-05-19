function [ res ] = attachUnit( value, unit )
%ATTACHUNIT Attaches unit to numerical value, returns COMSOL format string
%   Detailed explanation goes here
    res = [num2str(value),['[',unit,']']];
end

