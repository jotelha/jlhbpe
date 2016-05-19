function obj = update1dOperators(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% local definitions
%     fprintf('  Setting up local definitions...\n');
   
    surfaceVertexIndex = obj.m.selection('surfaceVertex1d').entities(0);
    zetaVertexIndex = obj.m.selection('zetaVertex1d').entities(0);
%     bulkExitVertexIndex = obj.m.selection('bulkExitVertex1d').entities(0);

    obj.m.cpl('projectReactionPlaneToSurface1d').selection.geom('geom1d', 0);
    obj.m.cpl('projectReactionPlaneToSurface1d').selection.named('zetaVertex1d');
    obj.m.cpl('projectReactionPlaneToSurface1d').label('projectReactionPlaneToSurface1d');
    obj.m.cpl('projectReactionPlaneToSurface1d').set('opname', 'projectReactionPlaneToSurface1d');
    obj.m.cpl('projectReactionPlaneToSurface1d').selection('dstvertex1').geom('geom1d', 0);
    obj.m.cpl('projectReactionPlaneToSurface1d').selection('dstvertex1').set(surfaceVertexIndex);
    obj.m.cpl('projectReactionPlaneToSurface1d').selection('srcvertex1').geom('geom1d', 0);
    obj.m.cpl('projectReactionPlaneToSurface1d').selection('srcvertex1').set(zetaVertexIndex);
 
    % integration operator on bpe surface
    obj.m.cpl('integrateSurface1d').selection.named('surfaceVertex1d');
    obj.m.cpl('integrateSurface1d').label('integrateSurface1d');
    obj.m.cpl('integrateSurface1d').set('opname', 'integrateSurface1d');

    % integrate exit towards bulk
    obj.m.cpl('integrateBulkExit1d').selection.named('bulkExitVertex1d');
    obj.m.cpl('integrateBulkExit1d').label('integrateBulkExit1d');
    obj.m.cpl('integrateBulkExit1d').set('opname', 'integrateBulkExit1d');

    % maximum on whole computational domain
%     obj.maximumOnDomain.set('opname', 'maximumOnDomain');
%     obj.maximumOnDomain.selection.all;
end