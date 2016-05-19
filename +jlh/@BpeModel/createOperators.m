function createOperators(obj)
    fprintf('  Creating local definitions...\n');
    %projectZetaPotential = obj.m.cpl.create('projectZetaPotential', 'LinearExtrusion', 'geom');
    obj.projectReactionPlaneToSurface = obj.m.cpl.create('projectReactionPlaneToSurface', 'LinearExtrusion', 'geom');

    % obj.projectSurfaceToBulkExit = obj.m.cpl.create('projectSurfaceToBulkExit', 'LinearExtrusion', 'geom');
    obj.integrateSurface = obj.m.cpl.create('integrateSurface', 'Integration', 'geom');
    obj.integrateBulkExit = obj.m.cpl.create('integrateBulkExit', 'Integration', 'geom');
    obj.maximumOnDomain = obj.m.cpl.create('maximumOnDomain', 'Maximum', 'geom');
end