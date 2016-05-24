function obj = updateMappedMesh(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    %% meshing
    fprintf('  Setting mesh parameters...\n');

     % get mesh refinement from parameters
%     obj.mesh1D();

    a0 = obj.intFirstDebyeLength(1);
    r = obj.intFirstDebyeLength(2)/obj.intFirstDebyeLength(1);
    N = ceil(log(1-obj.epsilon/a0*(1-r))/log(r));

    R = max(r^N, obj.intFirstDebyeLength(end)/obj.intFirstDebyeLength(1));

    a_lat = obj.intFirstDebyeLength(end);
    N_lat = ceil(obj.w_mesh / a_lat);
    
    N_lat_end = ceil((obj.w_mesh-obj.epsilon) / a_lat);

    hmax = obj.intRemaining(end);
    
    r_zeta = obj.intRemaining(1)/obj.intFirstDebyeLength(end);

    
    obj.m.mesh('mappedMesh').label('mappedMesh');
    
    obj.m.mesh('mappedMesh').feature('size').set('hmax', num2str(hmax) ); % is this the automatic refinement?
    obj.m.mesh('mappedMesh').feature('size').set('hmin', num2str(a0) ); 
    obj.m.mesh('mappedMesh').feature('size').set('hgrad', num2str(r_zeta) ); 

    % mesh refinement at BPE surface, whole surface including insulated
    % pieces
    
    obj.m.mesh('mappedMesh').feature('freeTriangularMesh').active(false); % only needed if more than 4 blocks
    obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').active(false);

    
    if obj.nMeshChops >= 4

        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').selection.geom('geom',1);
        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').selection.named('geom_unitCellPrototypeBpeSurface');
        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').label('edgeAtBpeSurface');

        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.geom('geom', 1); % what does this mean?
    %   obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('geom_unitCellPrototypeBpeSurface');
        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').label('bpeSurfaceMeshingProperties');
        obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').set('numelem', N_lat);


        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').selection.geom('geom',1);
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').selection.named('geom_unitCellPrototypeLowerLateralBoundary');
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').label('edgeAtFirstDebyeLengthLateralBoundary');

        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').selection.geom('geom', 1); % what does this mean?
    %   obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').selection.named('geom_unitCellPrototypeLowerLateralBoundary');
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').label('lateralBoundaryMeshingPropertiesAtFirstDebyeLength');
    %     obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('bpeSurfaceMeshingProperties').set('numelem', N_lat);
    %     obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('explicit', obj.distributionFirstDebyeLengthStr);
    %     obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('type', 'explicit');
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('type','predefined');
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('method', 'geometric');
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('elemcount', num2str(N));
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('elemratio', num2str(R));
        obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('reverse', 'on');


    %     growth = obj.distributionFirstDebyeLength(2:end)./ obj.distributionFirstDebyeLength(1:(end-1));
    %     obj.m.mesh('mappedMesh').feature('edgeAtFirstDebyeLengthLateralBoundary').feature('lateralBoundaryMeshingPropertiesAtFirstDebyeLength').set('hgrad', min(growth)); % maximum element growth coefficient

        % mesh size properties at zeta (or reaction) plane
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').selection.geom('geom',1);
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').selection.named('geom_unitCellPrototypeZetaPlane');
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').label('Mesh refinement at zeta plane');
    %     
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').selection.geom('geom', 1); % what does this mean?
    % %   obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').selection.named('geom_unitCellPrototypeZetaPlane');
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').label('Mesh refinement at zeta plane');
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('custom', 'on');
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hmaxactive', true);
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hmax', max(obj.intFirstDebyeLength)); % maximum element width
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgradactive', true);
    %     
    %     growth = obj.distributionFirstDebyeLength(2:end)./ obj.distributionFirstDebyeLength(1:(end-1));
    %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgrad', max(growth)); % maximum element growth coefficient

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
    % %   obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
    % %     obj.remainingBoundariesMeshingProperties.selection.remaining();
    %     obj.remainingBoundariesMeshingProperties.label('Mesh refinement for remaining boundaries');
    %     obj.remainingBoundariesMeshingProperties.set('custom', 'on');
    %     obj.remainingBoundariesMeshingProperties.set('hmaxactive', true);
    %     obj.remainingBoundariesMeshingProperties.set('hmax', max(obj.intRemaining)); % maximum element width
    % %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgradactive', true);
    % %     obj.m.mesh('mappedMesh').feature('zetaPlaneMeshingEdge').feature('zetaPlaneMeshingProperties').set('hgrad', max(obj.distributionFirstDebyeLength)); % maximum element growth coefficient

        % mesh in unit cell prototype 2d domain
        obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').active(true);
        obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').selection.geom('geom', 2);
        obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').selection.named('geom_unitCellPrototypeLowerDomain');
        obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').label('boundaryLayerMesh');

    %     obj.m.mesh('mappedMesh').feature('freeTriangularMesh').feature('domainMeshingProperties').selection.geom('geom', 2);
    %     obj.m.mesh('mappedMesh').feature('freeTriangularMesh').feature('domainMeshingProperties').selection.named('geom_unitCellPrototypeCumulative_dom');
    %     obj.m.mesh('mappedMesh').feature('freeTriangularMesh').feature('domainMeshingProperties').label('Free triangualar mesh maximum element size');
    %     obj.m.mesh('mappedMesh').feature('freeTriangularMesh').feature('domainMeshingProperties').set('custom', 'on');
    %     obj.m.mesh('mappedMesh').feature('freeTriangularMesh').feature('domainMeshingProperties').set('hmaxactive', true);
    %     obj.m.mesh('mappedMesh').feature('freeTriangularMesh').feature('domainMeshingProperties').set('hmax', hmax);
        obj.m.mesh('mappedMesh').feature('freeTriangularMesh').active(true);
        obj.m.mesh('mappedMesh').feature('freeTriangularMesh').selection.geom('geom', 2);
        obj.m.mesh('mappedMesh').feature('freeTriangularMesh').selection.named('geom_unitCellPrototypeUpperDomain');
        obj.m.mesh('mappedMesh').feature('freeTriangularMesh').label('freeTriangularMesh');
    end
    
    
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').selection.geom('geom',1);
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').selection.named('geom_eastwardThinningEdges');
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').label('eastwardThinningEdges');
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').selection.geom('geom', 1); % what does this mean?
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').selection.named('geom_eastwardThinningEdges');
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').label('eastwardThinningDistribution');
%     obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('type','predefined');
%     obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('method', 'geometric');
%     obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('elemcount', num2str(N));
%     obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('elemratio', num2str(R));
%     obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('reverse', 'off');
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('type','explicit');
    obj.m.mesh('mappedMesh').feature('eastwardThinningEdges').feature('eastwardThinningDistribution').set('explicit', obj.distributionFirstDebyeLengthStr);

    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').selection.geom('geom',1);
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').selection.named('geom_westwardThinningEdges');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').label('westwardThinningEdges');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').selection.geom('geom', 1); % what does this mean?
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').selection.named('geom_westwardThinningEdges');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').label('westwardThinningDistribution');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').set('type','predefined');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').set('method', 'geometric');
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').set('elemcount', num2str(N));
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').set('elemratio', num2str(R));
    obj.m.mesh('mappedMesh').feature('westwardThinningEdges').feature('westwardThinningDistribution').set('reverse', 'on');
    
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').selection.geom('geom',1);
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').selection.named('geom_lateralThinningEdges');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').label('lateralThinningEdges');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').selection.geom('geom', 1); % what does this mean?
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').selection.named('geom_lateralThinningEdges');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').label('lateralThinningDistribution');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').set('type','predefined');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').set('method', 'geometric');
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').set('elemcount', num2str(N));
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').set('elemratio', num2str(R));
    obj.m.mesh('mappedMesh').feature('lateralThinningEdges').feature('lateralThinningDistribution').set('reverse', 'off');
    
    % corrected edge vertex distribution at ends of bpe
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').selection.geom('geom',1);
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').selection.named('geom_bpeSurfaceAtBpeEnds');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').label('edgeAtBpeSurfaceAtBpeEnds');

    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').feature('bpeSurfaceAtBpeEndsMeshingProperties').selection.geom('geom', 1); % what does this mean?
%   obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').feature('bpeSurfaceAtBpeEndsMeshingProperties').selection.named('geom_bpeSurfaceAtBpeEnds');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').feature('bpeSurfaceAtBpeEndsMeshingProperties').label('bpeSurfaceAtBpeEndsMeshingProperties');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurfaceAtBpeEnds').feature('bpeSurfaceAtBpeEndsMeshingProperties').set('numelem', N_lat_end);

    obj.m.mesh('mappedMesh').feature('boundaryLayerMeshAtBpeEnds').selection.named('geom_mappedMeshDomainsAtEnds');
    obj.m.mesh('mappedMesh').feature('freeTriangularMeshAtBpeEnds').selection.named('geom_triangularMeshDomainsAtEnds');
%     obj.m.mesh('mappedMesh').feature('freeTriangularMeshAtDomainEnds').selection.named('geom_gapMeshingDomains');

%     obj.m.mesh('mappedMesh').run;
end