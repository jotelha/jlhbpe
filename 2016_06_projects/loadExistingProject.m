import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

%% restore BpeModel settings

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end


m = jlh.BpeModel;
m.projectName = tag;
m.projectPath = ['dat\',tag];
m.model_tag = tag;

fprintf('Calling %s...\n',caseStudyParameterFile);
run(caseStudyParameterFile);
m.prepareIdentifiers;

%% load model either from server or from file
tags = arrayfun(@(c) char(c),ModelUtil.tags,'UniformOutput',false); 
if any(strcmp(tags,m.model_tag))
    fprintf('Model %s exists on server, will be restored...\n',m.model_tag);
    m.m = ModelUtil.model(m.model_tag);
else
    mphFile = [m.projectPath,'\',m.projectName,'.mph'];
    fprintf('Model %s does not exist on server, loaded from file %s...',m.model_tag,mphFile);
    m.m = mphload(mphFile,m.model_tag);
end
model = m.m;
%% run logger
% also necessary to execute when loading

% kill logger, if already running
if exist('loggerPid','var') && loggerPid > 0
    ModelUtil.showProgress(false);
    fprintf('Trying to kill logger with PID %d...\n', loggerPid);
    system(sprintf('taskkill /pid %d /f',loggerPid),'-echo');
end

spawn = strrep('G:\scripts\launchers\spawn.bat','\','\\');
mTail = strrep('G:\scripts\mtail\mTail.exe','\','\\');

% server log file
logFile = [pwd(),'\',m.projectPath,'\comsol.log'];
% fclose(fopen(logFile, 'w'));
ModelUtil.showProgress(logFile);
logFileArg = strrep(logFile,'\','\\');

% system(sprintf('G:\\scripts\\mtail\\mTail.exe "%s" /start &',logFile));
cmd = prepTerm('spawn mTail "logFile" /start','spawn','mTail','logFile',spawn,mTail,logFileArg);
% spawn helper script yields logger pid
[status,cmdout] = system(cmd{1},'-echo');
loggerPid = str2double(cmdout);