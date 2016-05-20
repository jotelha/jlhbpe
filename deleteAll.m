fprintf('Removing project path %s\n',m.projectPath);
[status, message, messageid] = rmdir(m.projectPath,'s');
% system(sprintf('rmdir /S "%s"',m.projectPath),'-echo');
fprintf('%s.\n',message);
clearAll;
