function createChoppedMesh(obj)
    fprintf('  Creating mesh...\n');
    obj.standardMesh                        = obj.m.mesh.create('standardMesh', 'geom');    
    obj.meshingEdge                         = obj.standardMesh.create('meshingEdge', 'Edge'); % for 1d model
    obj.bpeSurfaceMeshingEdge               = obj.standardMesh.create('bpeSurfaceMeshingEdge', 'Edge');
    obj.bpeSurfaceMeshingProperties         = obj.bpeSurfaceMeshingEdge.create('bpeSurfaceMeshingProperties', 'Size');

    obj.zetaPlaneMeshingEdge                = obj.standardMesh.create('zetaPlaneMeshingEdge', 'Edge');
    obj.zetaPlaneMeshingProperties          = obj.zetaPlaneMeshingEdge.create('zetaPlaneMeshingProperties', 'Size');
    
%     obj.electrodeMeshingEdgeAtFirstDebyeLength          = obj.standardMesh.create('electrodeMeshingEdgeAtFirstDebyeLength', 'Edge');
    obj.lateralBoundaryMeshingPropertiesAtFirstDebyeLength   = obj.meshingEdge.create('lateralBoundaryMeshingPropertiesAtFirstDebyeLength', 'Distribution');
    obj.lateralBoundaryMeshingPropertiesAtRemainingDomain    = obj.meshingEdge.create('lateralBoundaryMeshingPropertiesAtRemainingDomain', 'Distribution');
   
%     obj.meshDistributionDDL                 = obj.meshingEdge.create('meshDistributionDDL', 'Distribution');
%     obj.meshDistributionRemaining           = obj.meshingEdge.create('meshDistributionRemaining', 'Distribution');
    %obj.electrodeMeshingPropertiesAtRemainingDomain     = obj.meshingEdge.create('electrodeMeshingPropertiesAtRemainingDomain', 'Distribution');
    %obj.remainingGeometryMeshingProperties = obj.meshingEdge.create('remainingGeometryMeshingProperties', 'Size');
   
    obj.remainingBoundariesMeshingEdge               = obj.standardMesh.create('remainingBoundariesMeshingEdge', 'Edge');
    obj.remainingBoundariesMeshingProperties         = obj.remainingBoundariesMeshingEdge.create('remainingBoundariesMeshingProperties', 'Size');

    obj.meshingDomain                       = obj.standardMesh.create('meshingDomain', 'FreeTri'); 
    obj.domainMeshingProperties             = obj.meshingDomain.create('domainMeshingProperties', 'Size');

   
    % a mesh to be filled "manually" here in MATLAB
    % obj.manualMesh = obj.m.mesh().create('manualMesh','geom');
    
    % a testing mesh
    obj.testingMesh                         = obj.m.mesh().create('manualMesh','geom');
    obj.testingMeshBpeSurfaceMeshingEdge    = obj.testingMesh.create('testingMeshBpeSurfaceMeshingEdge', 'Edge');
    obj.testingMeshBpeSurfaceMeshingProperties = obj.testingMeshBpeSurfaceMeshingEdge.create('testingMeshBpeSurfaceMeshingProperties', 'Size');

    obj.testingMeshTriangular               = obj.testingMesh.create('testingMeshTriangular', 'FreeTri'); 
    obj.testingMeshTriangularProperties     = obj.testingMeshTriangular.create('testingMeshTriangularProperties', 'Size');
end