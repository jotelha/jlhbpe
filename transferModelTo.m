function transferModelTo(m,server,port)
    import com.comsol.model.*
    import com.comsol.model.util.*
    currentlyOnline = true;
    if ~exist('server','var')
        server = 'localhost';
    end
    if ~exist('port','var')
        port = 2039;
    end
    
    try
    	tags = arrayfun(@(c) char(c),ModelUtil.tags,'UniformOutput',false); 
    catch
        tags = {};
        currentlyOnline = false;
    end
   
    % clean up currently running server
    % but do not remove model from server
    if currentlyOnline
        % save current model
        if ~isempty(m.model_tag)
            if any(strcmp(tags,m.model_tag))
                m.saveState;
            end
        end
        ModelUtil.showProgress(false);
        ModelUtil.disconnect();
    end
    
    % connect to other server
    mphstart(server,port)
    
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