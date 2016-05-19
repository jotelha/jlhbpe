%% all create commands and minimal settings
function obj = createModel(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    
    fprintf('  Creating model "%s" on server...\n',obj.model_tag);
    obj.m = ModelUtil.create(obj.model_tag);   % creates model on COMSOL server

    obj.createFunctions()
    
    % obj.comp_id = 'comp1';
    obj.m.modelNode.create(obj.comp_id); 

    obj.createGeometry();

    %% creating selections
    
    obj.createSelections();
    
    %% creating local definituins
    obj.createOperators()

    % variables
    obj.createVariables();

   %% create physics
    obj.createPhysics();
    %% creating mesh
    
    obj.createMesh();
    %% creating study
    %outsourced to "attachStudy"

    %% create probes
    obj.createProbes();
   %% create plots
   obj.create1dPlots();
   obj.create2dPlots();

      %% create views
   obj.standardView = obj.m.view.create('standardView', 'geom');
   obj.ddlView      = obj.m.view.create('ddlView', 'geom');
   
    %% create datasets
    obj.createDatasets();


end