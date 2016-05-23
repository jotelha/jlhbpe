if exist('loggerPid','var') && loggerPid > 0
    ModelUtil.showProgress(false);
    fprintf('Trying to kill logger with PID %d...\n', loggerPid);
    system(sprintf('taskkill /pid %d /f',loggerPid),'-echo');
end
modelTags = ModelUtil.tags;
for i = 1:numel(modelTags)
    ModelUtil.remove(modelTags(i));
end
ModelUtil.clear;
clear;
close all; % close all figures