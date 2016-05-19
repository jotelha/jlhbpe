function obj = clearProject(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    ModelUtil.remove(obj.model_tag);
    close all; % close all figures
    clear
end