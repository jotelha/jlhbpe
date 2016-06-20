function [titles,expressions,labels] = extractPlotInfo(plots)
	titles = plots.keys';
    values = plots.values';
    hasLabel = cellfun(@(c) iscell(c),values);
    labels = cell(numel(titles),1);
    labels(hasLabel) = cellfun(@(c) c{2}, values(hasLabel),'UniformOutput',false);
    values(hasLabel) = cellfun(@(c) c{1}, values(hasLabel),'UniformOutput',false);
    expressions = values;        
end