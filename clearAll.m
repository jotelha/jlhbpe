modelTags = ModelUtil.tags;
for i = numel(modelTags)
    ModelUtil.remove(modelTags(i));
end
ModelUtil.clear;
clear;
close all; % close all figures