function obj = updateChoppedMesh(obj)
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
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').selection.geom('geom',1);
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').selection.named('geom_unitCellPrototypeBpeSurface');
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').label('Mesh refinement at BPE surface');

    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').selection.geom('geom', 1); % what does this mean?
%   obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').selection.named('geom_unitCellPrototypeBpeSurface');
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').label('Mesh refinement at BPE surface');
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').set('custom', 'on');
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').set('hmaxactive', true);
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').set('hmax', min(obj.intFirstDebyeLength)); % maximum element width
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').set('hgradactive', true);
    
    growth = obj.distributionFirstDebyeLength(2:end)./ obj.distributionFirstDebyeLength(1:(end-1));
    obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').set('hgrad', min(growth)); % maximum element growth coefficient

    % mesh size properties at zeta (or reaction) plane
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').selection.geom('geom',1);
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').selection.named('geom_unitCellPrototypeZetaPlane');
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').label('Mesh refinement at zeta plane');
    
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').selection.geom('geom', 1); % what does this mean?
%   obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').selection.named('geom_unitCellPrototypeZetaPlane');
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').label('Mesh refinement at zeta plane');
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('custom', 'on');
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hmaxactive', true);
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hmax', max(obj.intFirstDebyeLength)); % maximum element width
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgradactive', true);
    
    growth = obj.distributionFirstDebyeLength(2:end)./ obj.distributionFirstDebyeLength(1:(end-1));
    obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgrad', max(growth)); % maximum element growth coefficient

    % mesh element distribution within diffuse layer
%     obj.electrodeMeshingEdgeAtFirstDebyeLength.selection.geom('geom',1);
%     obj.electrodeMeshingEdgeAtFirstDebyeLength.selection.named('electrodesAtFirstDebyeLength');
%     obj.electrodeMeshingEdgeAtFirstDebyeLength.label('Mesh refinement at electrodes at first debye length above bpe');
% 
%     obj.electrodeMeshingPropertiesAtFirstDebyeLength.selection.geom('geom', 1); % what does this mean?
%     obj.electrodeMeshingPropertiesAtFirstDebyeLength.selection.named('electrodesAtFirstDebyeLength');
%     obj.electrodeMeshingPropertiesAtFirstDebyeLength.label('Mesh refinement at electrodes at first debye length above bpe');
%     obj.electrodeMeshingPropertiesAtFirstDebyeLength.set('explicit', obj.distributionFirstDebyeLengthStr);
%     obj.electrodeMeshingPropertiesAtFirstDebyeLength.set('type', 'explicit');

%     obj.electrodeMeshingPropertiesAtRemainingDomain.active(false);
%     obj.electrodeMeshingPropertiesAtRemainingDomain.selection.geom('geom', 1); % what does this mean?
%     obj.electrodeMeshingPropertiesAtRemainingDomain.selection.named('electrodesAtRemainingDomain');
%     obj.electrodeMeshingPropertiesAtRemainingDomain.label('Mesh refinement at electrodes at remaining geometry');
%     obj.electrodeMeshingPropertiesAtRemainingDomain.set('explicit', obj.distributionRemainingStr);
%     obj.electrodeMeshingPropertiesAtRemainingDomain.set('type', 'explicit');
%     
    % meshing properties for remaining geometry
%     obj.remainingGeometryMeshingProperties.label('Mesh proterties for remaining geometry');
%     obj.remainingGeometryMeshingProperties.set('custom', 'on');
%     obj.remainingGeometryMeshingProperties.set('hmaxactive', true);
%     obj.remainingGeometryMeshingProperties.set('hmax', '1/1000'); % 1000 seems to be a good value for water
%     obj.remainingGeometryMeshingProperties.set('hgradactive', true);

    % mesh size properties for all other boundaries, only set maximum size
%     obj.remainingBoundariesMeshingEdge.selection.geom('geom',1);
%     obj.remainingBoundariesMeshingEdge.selection.remaining();
%     obj.remainingBoundariesMeshingEdge.label('Mesh refinement for remaining boundaries');
%     
% %     obj.remainingBoundariesMeshingProperties.selection.geom('geom', 1); % what does this mean?
% %   obj.m.mesh('standardMesh').feature('bpeSurfaceMeshingEdge').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
% %     obj.remainingBoundariesMeshingProperties.selection.remaining();
%     obj.remainingBoundariesMeshingProperties.label('Mesh refinement for remaining boundaries');
%     obj.remainingBoundariesMeshingProperties.set('custom', 'on');
%     obj.remainingBoundariesMeshingProperties.set('hmaxactive', true);
%     obj.remainingBoundariesMeshingProperties.set('hmax', max(obj.intRemaining)); % maximum element width
% %     obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgradactive', true);
% %     obj.m.mesh('standardMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgrad', max(obj.distributionFirstDebyeLength)); % maximum element growth coefficient

    % mesh in unit cell prototype 2d domain
    obj.m.mesh('standardMesh').feature('meshingDomain').selection.geom('geom', 2);
    obj.m.mesh('standardMesh').feature('meshingDomain').selection.named('geom_unitCellPrototypeCumulative_dom');
    obj.m.mesh('standardMesh').feature('meshingDomain').label('Free triangualar mesh');

    obj.m.mesh('standardMesh').feature('meshingDomain').feature('domainMeshingProperties').selection.geom('geom', 2);
    obj.m.mesh('standardMesh').feature('meshingDomain').feature('domainMeshingProperties').selection.named('geom_unitCellPrototypeCumulative_dom');
    obj.m.mesh('standardMesh').feature('meshingDomain').feature('domainMeshingProperties').label('Free triangualar mesh maximum element size');
    obj.m.mesh('standardMesh').feature('meshingDomain').feature('domainMeshingProperties').set('custom', 'on');
    obj.m.mesh('standardMesh').feature('meshingDomain').feature('domainMeshingProperties').set('hmaxactive', true);
    obj.m.mesh('standardMesh').feature('meshingDomain').feature('domainMeshingProperties').set('hmax', max(obj.intRemaining));
    
    % testing mesh
    obj.testingMesh.label('testingMesh');
    %obj.standardMesh.feature('size').set('hauto', 2); % is this the automatic refinement?
   
    % mesh refinement at BPE surface, whole surface including insulated
    % pieces
    obj.testingMeshBpeSurfaceMeshingEdge.selection.geom('geom',1);
    obj.testingMeshBpeSurfaceMeshingEdge.selection.named('geom_unitCellPrototypeBpeSurface');
    obj.testingMeshBpeSurfaceMeshingEdge.label('Mesh refinement at BPE surface');

    obj.testingMeshBpeSurfaceMeshingProperties.selection.geom('geom', 1); % what does this mean?
    obj.testingMeshBpeSurfaceMeshingProperties.selection.named('geom_unitCellPrototypeBpeSurface');
    obj.testingMeshBpeSurfaceMeshingProperties.label('Mesh refinement at BPE surface');
    obj.testingMeshBpeSurfaceMeshingProperties.set('custom', 'on');
    obj.testingMeshBpeSurfaceMeshingProperties.set('hmaxactive', true);
    obj.testingMeshBpeSurfaceMeshingProperties.set('hmax', 'epsilon/10'); % maximum element width
    obj.testingMeshBpeSurfaceMeshingProperties.set('hgradactive', true);
    obj.testingMeshBpeSurfaceMeshingProperties.set('hgrad', '1.01'); % maximum element growth coefficient

    % mesh in remaining 2d domain
    obj.testingMeshTriangular.selection.geom('geom', 2);
    obj.testingMeshTriangular.selection.named('geom_unitCellPrototypeCumulative_dom');
    obj.testingMeshTriangular.label('Free triangualar mesh for testing purposes');
    
    obj.testingMeshTriangularProperties.selection.geom('geom', 2);
    obj.testingMeshTriangularProperties.selection.named('geom_unitCellPrototypeCumulative_dom');
    obj.testingMeshTriangularProperties.label('Free triangualar mesh for testing purposes maximum element size');
    obj.testingMeshTriangularProperties.set('custom', 'on');
    obj.testingMeshTriangularProperties.set('hmaxactive', true);
    obj.testingMeshTriangularProperties.set('hmax', '1/4');
    % obj.standardMesh.run;
end