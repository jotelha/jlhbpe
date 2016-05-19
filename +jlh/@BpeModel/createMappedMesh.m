function createMappedMesh(obj)
% mapped mesh
    obj.m.mesh().create('mappedMesh','geom');
    obj.m.mesh('mappedMesh').create('edgeAtFirstDebyeLengthLateralBoundary','Edge');
    obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').create('lateralBoundaryMeshingPropertiesAtFirstDebyeLength','Distribution');

    obj.m.mesh('mappedMesh').create('edgeAtBpeSurface','Edge');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').create('bpeSurfaceMeshingProperties','Distribution');

    obj.m.mesh('mappedMesh').create('westwardThinningEdges','Edge');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').create('westwardThinningDistribution','Distribution');
    
    obj.m.mesh('mappedMesh').create('eastwardThinningEdges','Edge');
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').create('eastwardThinningDistribution','Distribution');
    
    obj.m.mesh('mappedMesh').create('lateralThinningEdges','Edge');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').create('lateralThinningDistribution','Distribution');
    
    obj.m.mesh('mappedMesh').create('edgeAtBpeSurfaceAtBpeEnds','Edge');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').create('bpeSurfaceAtBpeEndsMeshingProperties','Distribution');

    obj.m.mesh('mappedMesh').create('boundaryLayerMesh','Map');
    obj.m.mesh('mappedMesh').create('freeTriangularMesh','FreeTri');
    
    obj.m.mesh('mappedMesh').create('boundaryLayerMeshAtBpeEnds','Map');
    obj.m.mesh('mappedMesh').create('freeTriangularMeshAtBpeEnds','FreeTri');
%     obj.m.mesh('mappedMesh').create('freeTriangularMeshAtDomainEnds','FreeTri');

end