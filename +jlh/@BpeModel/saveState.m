%% saves the current model and important state variables 
function obj = saveState(obj,name)
    if isempty(obj.projectName)
        obj.newProject(name);
    end

%     projectPath = ['dat/',projectName];

    fprintf('Saving project %s to directory %s...\n',obj.projectName,obj.projectPath);

    mphFile = [obj.projectPath,'/',obj.projectName,'.mph'];
    matFile = [obj.projectPath,'/',obj.projectName,'.mat'];
%     datFile = [projectPath,'/',projectName,'_dat.mat'];

    if( ~exist(obj.projectPath,'dir') )
        fprintf('Creating project directory %s...\n',obj.projectPath);
        mkdir(obj.projectPath);
    end

    fprintf('Copying plots and parameter files...\n');
%     copyfile('img/*.png', projectPath);
%     copyfile('setGeneralParameters.m', projectPath);
%     copyfile('setReactionParameters.m', projectPath);

    fprintf('Saving MATLAB class to %s...\n',matFile);
    save(matFile,'obj');

    %% saves data of all iteration steps
%     fprintf('Extracting solution data from COMSOL model and saving to %s...\n',datFile);
%     extractAllData
%     save(datFile,'p','phi','c');

    %% save COMSOL model
    fprintf('Saving COMSOL model to %s...\n',mphFile);
    mphsave(obj.m,mphFile);
end
