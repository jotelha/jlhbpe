function saveMapAsTxt(map,txt)
    %% save file
    key = map.keys';
    value = map.values';
    hasDesc = cellfun(@(c) iscell(c),value);
    comment = cell(numel(key),1);
    comment(hasDesc) = cellfun(@(c) c{2}, value(hasDesc),'UniformOutput',false);
    value(hasDesc) = cellfun(@(c) c{1}, value(hasDesc),'UniformOutput',false);

    if all(cellfun(@(c) isempty(c),comment))
        T = table(key,value);
    else
        T = table(key,value,comment);
    end
    writetable(T,txt,'WriteRowNames',false,'WriteVariableNames',false,'Delimiter',' ');
end