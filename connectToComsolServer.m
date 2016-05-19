function connectToComsolServer(server,port)
    if ~exist('port','var')
        port = 2036;
    end
    if ~exist('server','var')
        server = '166.111.52.27';
    end
    mphstart(server,port)
    import com.comsol.model.*
    import com.comsol.model.util.*
end