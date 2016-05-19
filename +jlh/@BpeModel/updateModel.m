% 
% all model settings except create commands
%
function obj = updateModel(obj)
%     import com.comsol.model.*
%     import com.comsol.model.util.*
%     import jlh.hf.*
%     import jlh.*
    
    fprintf('  Modify settings of model "%s"...\n',obj.model_tag);

    obj.updateParameters();
    obj.updateFunctions();
    obj.updateGeometry();
    obj.updateSelections();
    obj.updateOperators();
    obj.updatePhysics();
    obj.updateMesh();
end