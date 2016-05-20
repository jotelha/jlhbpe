function switchServerTo(server,port)
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
    if currentlyOnline
        ModelUtil.showProgress(false);
        ModelUtil.disconnect();
    end
    
    % connect to other server
    mphstart(server,port)
end