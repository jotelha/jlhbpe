function obj = updateOperatorsForChoppedGeometry(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% local definitions
    fprintf('  Setting up local definitions...\n');
    % % operators
    % projectZetaPotential.selection.geom('geom', 0);
    % projectZetaPotential.selection.named('zetaNode');
    % projectZetaPotential.label('projectZetaPotential');
    % projectZetaPotential.set('opname', 'projectZetaPotential');
    % projectZetaPotential.selection('dstvertex1').geom('geom', 0);
    % projectZetaPotential.selection('dstvertex1').set(1);
    % %projectZetaPotential.selection('dstvertex1').named('surfaceNode');
    % %projectZetaPotential.selection('dstvertex2').geom('geom', 0);
    % %projectZetaPotential.selection('dstvertex3').geom('geom', 0);
    % projectZetaPotential.selection('srcvertex1').geom('geom', 0);
    % projectZetaPotential.selection('srcvertex1').set(2);
    % %projectZetaPotential.selection('srcvertex1').named('zetaNode');
    % %projectZetaPotential.selection('srcvertex2').geom('geom', 0);
    % %projectZetaPotential.selection('srcvertex3').geom('geom', 0);

    % operators
    
    leftBoundaryOfSurfaceVertexIndex = obj.leftBoundaryOfSurface.entities(0);
    rightBoundaryOfSurfaceVertexIndex = obj.rightBoundaryOfSurface.entities(0);
    leftBoundaryOfZetaPlaneVertexIndex = obj.leftBoundaryOfZetaPlane.entities(0);
    rightBoundaryOfZetaPlaneVertexIndex = obj.rightBoundaryOfZetaPlane.entities(0);

    obj.projectReactionPlaneToSurface.selection.geom('geom', 0);
    obj.projectReactionPlaneToSurface.selection.named('geom_zetaPlane');
    obj.projectReactionPlaneToSurface.label('projectReactionPlaneToSurface');
    obj.projectReactionPlaneToSurface.set('opname', 'projectReactionPlaneToSurface');
%     obj.projectReactionPlaneToSurface.selection('dstvertex1').geom('geom', 0);
    obj.projectReactionPlaneToSurface.selection('dstvertex1').set(leftBoundaryOfSurfaceVertexIndex(1));
    obj.projectReactionPlaneToSurface.selection('dstvertex2').set(rightBoundaryOfSurfaceVertexIndex(1));
%     obj.projectReactionPlaneToSurface.selection('srcvertex1').geom('geom', 0);
    obj.projectReactionPlaneToSurface.selection('srcvertex1').set(leftBoundaryOfZetaPlaneVertexIndex(1));
    obj.projectReactionPlaneToSurface.selection('srcvertex2').set(rightBoundaryOfZetaPlaneVertexIndex(1));

    
    % project currents at surface to bulk exit
%     obj.projectSurfaceToBulkExit.selection.geom('geom', 0);
%     obj.projectSurfaceToBulkExit.selection.named('bulkBoundary');
%     obj.projectSurfaceToBulkExit.label('projectSurfaceToBulkExit');
%     obj.projectSurfaceToBulkExit.set('opname', 'projectSurfaceToBulkExit');
%     obj.projectSurfaceToBulkExit.selection('dstvertex1').geom('geom', 0);
%     obj.projectSurfaceToBulkExit.selection('dstvertex1').set(3);
%     obj.projectSurfaceToBulkExit.selection('srcvertex1').geom('geom', 0);
%     obj.projectSurfaceToBulkExit.selection('srcvertex1').set(1);
%     
    % integration operator on bpe surface
    obj.integrateSurface.selection.named('geom_bpeSurface');
    obj.integrateSurface.label('integrateSurface');
    obj.integrateSurface.set('opname', 'integrateSurface');

    % integrate exit towards bulk
    obj.integrateBulkExit.selection.named('geom_bulkBoundary');
    obj.integrateBulkExit.label('integrateBulkBoundary');
    obj.integrateBulkExit.set('opname', 'integrateBulkBoundary');

    % maximum on whole computational domain
    obj.maximumOnDomain.set('opname', 'maximumOnDomain');
    obj.maximumOnDomain.selection.all;
end