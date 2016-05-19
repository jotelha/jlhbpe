function create1dFluxOperators(obj)
    fprintf('  Creating local definitions...\n');

    % obj.projectSurfaceToBulkExit = obj.m.cpl.create('projectSurfaceToBulkExit', 'LinearExtrusion', 'geom');
    obj.integrateWE = obj.m.cpl.create('integrateWE', 'Integration', 'geom');
    obj.integrateCE = obj.m.cpl.create('integrateCE', 'Integration', 'geom');
    obj.integrateDomain = obj.m.cpl.create('integrateDomain', 'Integration', 'geom');

    obj.maximumOnDomain = obj.m.cpl.create('maximumOnDomain', 'Maximum', 'geom');
    
end