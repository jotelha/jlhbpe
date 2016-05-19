function obj = updateMeshMapped(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    %% meshing
    fprintf('  Setting mesh parameters...\n');

     % get mesh refinement from parameters
    obj.mesh1D();
    
    obj.standardMesh.label('standardMesh');
    %obj.standardMesh.feature('size').set('hauto', 2); % is this the automatic refinement?
   
    % mesh refinement at BPE surface, whole surface including insulated
    % pieces
    obj.bpeSurfaceMeshingEdge.selection.geom('geom',1);
    obj.bpeSurfaceMeshingEdge.selection.named('entireSurface');
    obj.bpeSurfaceMeshingEdge.label('Mesh refinement at BPE surface');

    obj.bpeSurfaceMeshingProperties.selection.geom('geom', 1); % what does this mean?
%   obj.bpeSurfaceMeshingProperties.selection.named('bpeSurface');
    obj.bpeSurfaceMeshingProperties.selection.named('entireSurface');
    obj.bpeSurfaceMeshingProperties.label('Mesh refinement at BPE surface');
    obj.bpeSurfaceMeshingProperties.set('custom', 'on');
    obj.bpeSurfaceMeshingProperties.set('hmaxactive', true);
    obj.bpeSurfaceMeshingProperties.set('hmax', min(obj.intFirstDebyeLength)); % maximum element width
    obj.bpeSurfaceMeshingProperties.set('hgradactive', true);
    
    growth = obj.distributionFirstDebyeLength(2:end)./ obj.distributionFirstDebyeLength(1:(end-1));
    obj.bpeSurfaceMeshingProperties.set('hgrad', min(growth)); % maximum element growth coefficient

    % mesh size properties at zeta (or reaction) plane
    obj.zetaPlaneMeshingEdge.selection.geom('geom',1);
    obj.zetaPlaneMeshingEdge.selection.named('zetaPlane');
    obj.zetaPlaneMeshingEdge.label('Mesh refinement at zeta plane');
    
    obj.zetaPlaneMeshingProperties.selection.geom('geom', 1); % what does this mean?
%   obj.bpeSurfaceMeshingProperties.selection.named('bpeSurface');
    obj.zetaPlaneMeshingProperties.selection.named('zetaPlane');
    obj.zetaPlaneMeshingProperties.label('Mesh refinement at zeta plane');
    obj.zetaPlaneMeshingProperties.set('custom', 'on');
    obj.zetaPlaneMeshingProperties.set('hmaxactive', true);
    obj.zetaPlaneMeshingProperties.set('hmax', max(obj.intFirstDebyeLength)); % maximum element width
    obj.zetaPlaneMeshingProperties.set('hgradactive', true);
    
    growth = obj.distributionFirstDebyeLength(2:end)./ obj.distributionFirstDebyeLength(1:(end-1));
    obj.zetaPlaneMeshingProperties.set('hgrad', max(growth)); % maximum element growth coefficient

    % mesh element distribution within diffuse layer
   
    obj.meshingEdge.selection.geom('geom',1');
    obj.meshingEdge.selection.named('lateralBoundary');

    
%     obj.electrodeMeshingEdgeAtFirstDebyeLength.selection.geom('geom',1);
%     obj.electrodeMeshingEdgeAtFirstDebyeLength.selection.named('lateralBoundaryAtFirstDebyeLength');
%     obj.electrodeMeshingEdgeAtFirstDebyeLength.label('Mesh refinement at electrodes at first debye length above bpe');

    obj.lateralBoundaryMeshingPropertiesAtFirstDebyeLength.selection.geom('geom', 1); % what does this mean?
    obj.lateralBoundaryMeshingPropertiesAtFirstDebyeLength.selection.named('lateralBoundaryAtFirstDebyeLength');
    obj.lateralBoundaryMeshingPropertiesAtFirstDebyeLength.label('Mesh refinement at electrodes beyond first debye length above bpe');
    obj.lateralBoundaryMeshingPropertiesAtFirstDebyeLength.set('explicit', obj.distributionFirstDebyeLengthStr);
    obj.lateralBoundaryMeshingPropertiesAtFirstDebyeLength.set('type', 'explicit');

%     obj.lateralBoundaryMeshingPropertiesAtRemainingDomain.active(false);
    obj.lateralBoundaryMeshingPropertiesAtRemainingDomain.selection.geom('geom', 1); % what does this mean?
    obj.lateralBoundaryMeshingPropertiesAtRemainingDomain.selection.named('lateralBoundaryAtRemainingDomain');
    obj.lateralBoundaryMeshingPropertiesAtRemainingDomain.label('Mesh refinement at electrodes at remaining geometry');
    obj.lateralBoundaryMeshingPropertiesAtRemainingDomain.set('explicit', obj.distributionRemainingStr);
    obj.lateralBoundaryMeshingPropertiesAtRemainingDomain.set('type', 'explicit');
    
    % meshing properties for remaining geometry
%     obj.remainingGeometryMeshingProperties.label('Mesh proterties for remaining geometry');
%     obj.remainingGeometryMeshingProperties.set('custom', 'on');
%     obj.remainingGeometryMeshingProperties.set('hmaxactive', true);
%     obj.remainingGeometryMeshingProperties.set('hmax', '1/1000'); % 1000 seems to be a good value for water
%     obj.remainingGeometryMeshingProperties.set('hgradactive', true);

    % mesh size properties for all other boundaries, only set maximum size
    obj.remainingBoundariesMeshingEdge.selection.geom('geom',1);
    obj.remainingBoundariesMeshingEdge.selection.remaining();
    obj.remainingBoundariesMeshingEdge.label('Mesh refinement for remaining boundaries');
    
%     obj.remainingBoundariesMeshingProperties.selection.geom('geom', 1); % what does this mean?
%   obj.bpeSurfaceMeshingProperties.selection.named('bpeSurface');
%     obj.remainingBoundariesMeshingProperties.selection.remaining();
    obj.remainingBoundariesMeshingProperties.label('Mesh refinement for remaining boundaries');
    obj.remainingBoundariesMeshingProperties.set('custom', 'on');
    obj.remainingBoundariesMeshingProperties.set('hmaxactive', true);
    obj.remainingBoundariesMeshingProperties.set('hmax', max(obj.intRemaining)); % maximum element width
%     obj.zetaPlaneMeshingProperties.set('hgradactive', true);
%     obj.zetaPlaneMeshingProperties.set('hgrad', max(obj.distributionFirstDebyeLength)); % maximum element growth coefficient

    % mesh in remaining 2d domain
    obj.meshingDomain.selection.geom('geom', 2);
    obj.meshingDomain.selection.all;
    obj.meshingDomain.label('Free triangualar mesh');

    
    obj.domainMeshingProperties.selection.geom('geom', 2);
    obj.domainMeshingProperties.selection.all;
    obj.domainMeshingProperties.label('Free triangualar mesh maximum element size');
    obj.domainMeshingProperties.set('custom', 'on');
    obj.domainMeshingProperties.set('hmaxactive', true);
    obj.domainMeshingProperties.set('hmax', max(obj.intRemaining));
    
    % testing mesh
    obj.testingMesh.label('testingMesh');
    %obj.standardMesh.feature('size').set('hauto', 2); % is this the automatic refinement?
   
    % mesh refinement at BPE surface, whole surface including insulated
    % pieces
    obj.testingMeshBpeSurfaceMeshingEdge.selection.geom('geom',1);
    obj.testingMeshBpeSurfaceMeshingEdge.selection.named('entireSurface');
    obj.testingMeshBpeSurfaceMeshingEdge.label('Mesh refinement at BPE surface');

    obj.testingMeshBpeSurfaceMeshingProperties.selection.geom('geom', 1); % what does this mean?
    obj.testingMeshBpeSurfaceMeshingProperties.selection.named('entireSurface');
    obj.testingMeshBpeSurfaceMeshingProperties.label('Mesh refinement at BPE surface');
    obj.testingMeshBpeSurfaceMeshingProperties.set('custom', 'on');
    obj.testingMeshBpeSurfaceMeshingProperties.set('hmaxactive', true);
    obj.testingMeshBpeSurfaceMeshingProperties.set('hmax', 'epsilon/10'); % maximum element width
    obj.testingMeshBpeSurfaceMeshingProperties.set('hgradactive', true);
    obj.testingMeshBpeSurfaceMeshingProperties.set('hgrad', '1.01'); % maximum element growth coefficient

    % mesh in remaining 2d domain
    obj.testingMeshTriangular.selection.geom('geom', 2);
    obj.testingMeshTriangular.selection.all;
    obj.testingMeshTriangular.label('Free triangualar mesh for testing purposes');
    
    obj.testingMeshTriangularProperties.selection.geom('geom', 2);
    obj.testingMeshTriangularProperties.selection.all;
    obj.testingMeshTriangularProperties.label('Free triangualar mesh for testing purposes maximum element size');
    obj.testingMeshTriangularProperties.set('custom', 'on');
    obj.testingMeshTriangularProperties.set('hmaxactive', true);
    obj.testingMeshTriangularProperties.set('hmax', '1/4');
    % obj.standardMesh.run;
end