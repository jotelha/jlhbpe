% assigns unique model tag and project name with date and time
function obj = newProject(obj,name)
    obj.projectName = [datestr(clock,'yyyy_mm_dd_HH_MM_SS_'),name];
    % model_tag = 'equilibrateAndSweepSurfacePotential';
    obj.model_tag = obj.projectName;
    % tempPath = 'tmp';
    obj.projectPath = ['dat\',obj.projectName];
    if( ~exist(obj.projectPath,'dir') )
        fprintf('Creating project directory %s...\n',obj.projectPath);
        mkdir(obj.projectPath);
    end
end