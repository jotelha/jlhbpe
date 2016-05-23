modelTags = ModelUtil.tags;
for i = 1:numel(modelTags)
    ModelUtil.remove(modelTags(i));
end
ModelUtil.clear;
clear;
close all; % close all figures