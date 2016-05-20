fprintf('Removing project path %s.\n',m.projectPath);
[status, message, messageid] = rmdir(m.projectPath,'s');
fprintf('%s.\n',message);
clearAll;
