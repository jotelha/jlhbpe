function obj = updateMappedMeshForSimpleGeometry(obj)
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

%     a_lat = obj.intFirstDebyeLength(end);
%     N_lat = ceil(obj.w_mesh / a_lat);
    a_lat = obj.intFirstDebyeLength(end);
    N_lat = ceil(obj.w / a_lat);
    
%     N_lat_end = ceil((obj.w_mesh-obj.epsilon) / a_lat);

    hmax = obj.intRemaining(end);
    
    r_zeta = obj.intRemaining(1)/obj.intFirstDebyeLength(end);

    
    obj.m.mesh('mappedMesh').label('mappedMesh');
    
    obj.m.mesh('mappedMesh').feature('size').set('hmax', num2str(hmax) ); % is this the automatic refinement?
    obj.m.mesh('mappedMesh').feature('size').set('hmin', num2str(a0) ); 
    obj.m.mesh('mappedMesh').feature('size').set('hgrad', num2str(r_zeta) ); 
    
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
    
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').selection.geom('geom',1);
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').selection.named('geom_entireSurface');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').label('edgeAtBpeSurface');

    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.geom('geom', 1); % what does this mean?
%   obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('bpeSurface');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').selection.named('geom_entireSurface');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').label('bpeSurfaceMeshingProperties');
    obj.m.mesh('mappedMesh').feature('edgeAtBpeSurface').feature('bpeSurfaceMeshingProperties').set('numelem', N_lat);
    
    obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').selection.geom('geom', 2);
    obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').selection.named('geom_ddl_dom');
    obj.m.mesh('mappedMesh').feature('boundaryLayerMesh').label('boundaryLayerMesh');

    obj.m.mesh('mappedMesh').feature('freeTriangularMesh').selection.geom('geom', 2);
    obj.m.mesh('mappedMesh').feature('freeTriangularMesh').selection.named('geom_space_dom');
    obj.m.mesh('mappedMesh').feature('freeTriangularMesh').label('freeTriangularMesh');

%     obj.m.mesh('mappedMesh').run;
end