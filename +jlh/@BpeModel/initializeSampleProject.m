function obj = initializeSampleProject(obj)
    import jlh.hf.*
    close all;
    
    obj.newProject('sample_project');
        
%     obj.setKrawiecCaseGeneralParameters();
%     obj.setKrawiecCaseReactionParameters();
    obj.setSampleCaseParameters();

    obj.prepareIdentifiers();
    obj.createModel();
    
    obj.updateModel();
    obj.addParametricSweepStudy('chargeDensityRampFactor',0);
end