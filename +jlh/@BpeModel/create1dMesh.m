function create1dMesh(obj)
    fprintf('  Creating mesh...\n');
    obj.standardMesh                        = obj.m.mesh.create('standardMesh', 'geom');    
    obj.meshingEdge                         = obj.standardMesh.create('meshingEdge', 'Edge'); % for 1d model
    obj.meshDistributionDDL                 = obj.meshingEdge.create('meshDistributionDDL', 'Distribution');
    obj.meshDistributionRemaining           = obj.meshingEdge.create('meshDistributionRemaining', 'Distribution');
end