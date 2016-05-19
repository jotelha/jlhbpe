function update1dModel(obj)
    obj.updateParameters();
    obj.updateFunctions();
    
    obj.update1dGeometry();
    obj.update1dSelections();
    obj.update1dOperators();
    obj.update1dPhysics();
    obj.update1dMesh();
end
    
