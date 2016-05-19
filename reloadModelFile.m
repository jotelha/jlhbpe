function reloadModelFile(m)
    import com.comsol.model.*
    import com.comsol.model.util.*
   
    % restore logging
    logFile = [pwd(),'\',m.projectPath,'\comsol.log'];
    % fclose(fopen(logFile, 'w'));
    ModelUtil.showProgress(logFile);
        
    if isempty(m.model_tag)
        m.model_tag = ModelUtil.uniquetag();
    else
        % if this model exists on server,
        % replace by currently transferred model
        tags = arrayfun(@(c) char(c),ModelUtil.tags,'UniformOutput',false); 
        if any(strcmp(tags,m.model_tag))
            ModelUtil.remove(m.model_tag);
        end
    end

    % load model to server  
    mphFile = [m.projectPath,'/',m.projectName,'.mph'];
    m.m = mphload(mphFile,m.model_tag);
end