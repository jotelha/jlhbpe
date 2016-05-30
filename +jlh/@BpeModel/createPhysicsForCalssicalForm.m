function obj = createPhysicsForClassicalForm(obj)
     %% creating physics
    obj.m.physics.create('NernstPlanckEquation','DilutedSpecies','geom',obj.c_id);
    
    obj.m.physics('NernstPlanckEquation').feature.create('zetaPlaneContinuity', 'Continuity', 1);

    obj.m.physics('NernstPlanckEquation').feature.create('FluxAtSurface', 'Fluxes', 1);
    obj.m.physics('NernstPlanckEquation').feature.create('BulkConcentrationsBC', 'Concentration', 1);

    obj.m.physics.create('PoissonEquation', 'Electrostatics', 'geom');

    obj.PoissonEquation.create('bpeChargeDensityBC', 'SurfaceChargeDensity', 1);     
%     obj.PoissonEquation.create('bpePotentialBC', 'DirichletBoundary', 1);
    
%     obj.bulkPotentialBC = obj.PoissonEquation.create('bulkPotentialBC', 'ElectricPotential', 1);
    obj.PoissonEquation.create('electrodePotentialBC', 'ElectricPotential', 1);  

    obj.zeroSurfaceCurrent = obj.m.physics.create('zeroSurfaceCurrent', 'GlobalEquations', 'geom');
    obj.zeroSurfaceCurrentEquation = obj.zeroSurfaceCurrent.create('zeroSurfaceCurrentEquation', 'GlobalEquations');
end