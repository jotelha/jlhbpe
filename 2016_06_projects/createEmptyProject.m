import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

if ~exist('files','var')
    files = containers.Map;
end

%% create new project

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end

m = jlh.BpeModel;
fprintf('Calling %s...\n',caseStudyParameterFile);
run(caseStudyParameterFile);
m.prepareIdentifiers;

if ~exist('caseStudyTitle','var')
    caseStudyTitle = 'emptyProject';
end
if ~exist('caseStudyTitleSuffix','var')
    caseStudyTitleSuffix = '';
end

m.newProject([caseStudyTitle,caseStudyTitleSuffix]);

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


%% remove other models from server, create COMSOL empty model
loadedModels    = ModelUtil.tags;
isLoaded        = arrayfun( @(s) strcmp(m.model_tag,s),loadedModels);
if any(isLoaded)
    ModelUtil.remove(m.model_tag);
end
m.m = ModelUtil.create(m.model_tag);   % creates model on COMSOL server