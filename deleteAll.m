if exist('loggerPid','var') && loggerPid > 0
    ModelUtil.showProgress(false);
    fprintf('Trying to kill logger with PID %d...\n', loggerPid);
    system(sprintf('taskkill /pid %d /f',loggerPid),'-echo');
    clear loggedPid;
end
fprintf('Removing project path %s\n',m.projectPath);
% [status, message, messageid] = rmdir(m.projectPath,'s');
system(sprintf('rmdir /s /q "%s"',m.projectPath),'-echo');
% fprintf('%s.\n',message);
clearAll;
