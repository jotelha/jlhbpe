function writeParameterTextFile(parameters,file)
    parameter = parameters.keys';
    value = parameters.values';
    hasDesc = cellfun(@(c) iscell(c) && numel(c) > 1,value);
    comment = cell(numel(parameter),1);
    comment(hasDesc) = cellfun(@(c) c{2}, value(hasDesc),'UniformOutput',false);
    value(hasDesc) = cellfun(@(c) c{1}, value(hasDesc),'UniformOutput',false);

    T = table(parameter,value,comment);
    writetable(T,file,'WriteRowNames',false,'WriteVariableNames',false,'Delimiter',' ');
end