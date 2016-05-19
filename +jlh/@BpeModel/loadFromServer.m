function obj = loadFromServer(obj,tag)
    import com.comsol.model.m.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    obj.m = ModelUtil.model(tag);
end