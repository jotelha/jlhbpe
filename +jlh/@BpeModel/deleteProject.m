function obj = deleteProject(obj)
    obj.clearProject;
    if exist(obj.projectPath,'dir');
        rmdir(obj.projectPath,'s');
    end
    obj = jlh.BpeModel;
end