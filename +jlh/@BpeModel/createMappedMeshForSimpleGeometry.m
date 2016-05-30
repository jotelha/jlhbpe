function createMappedMeshForSimpleGeometry(obj)
% mapped mesh
    obj.m.mesh().create('mappedMesh','geom');
%     obj.m.mesh('mappedMeshForSimpleGeometry').create('edgeAtFirstDebyeLengthLateralBoundary','Edge');
%     obj.m.mesh('mappedMeshForSimpleGeometry').feature('edgeAtFirstDebyeLengthLateralBoundary').create('lateralBoundaryMeshingPropertiesAtFirstDebyeLength','Distribution');

    obj.m.mesh('mappedMesh').create('edgeAtBpeSurface','Edge');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').create('bpeSurfaceMeshingProperties','Distribution');

    obj.m.mesh('mappedMesh').create('lateralThinningEdges','Edge');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').create('lateralThinningDistribution','Distribution');
    
%     obj.m.mesh('mappedMesh').create('edgeAtBpeSurfaceAtBpeEnds','Edge');
%     obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').create('bpeSurfaceAtBpeEndsMeshingProperties','Distribution');

    obj.m.mesh('mappedMesh').create('boundaryLayerMesh','Map');
    obj.m.mesh('mappedMesh').create('freeTriangularMesh','FreeTri');
    
%     obj.m.mesh('mappedMesh').create('boundaryLayerMeshAtBpeEnds','Map');
%     obj.m.mesh('mappedMesh').create('freeTriangularMeshAtBpeEnds','FreeTri');
%     obj.m.mesh('mappedMesh').create('freeTriangularMeshAtDomainEnds','FreeTri');

end