function obj = createPhysicsForClassicalForm(obj)
     %% creating physics
    obj.m.physics.create('NernstPlanckEquation','DilutedSpecies','geom',obj.c_id);
    
    obj.m.physics('NernstPlanckEquation').feature.create('FluxAtSurface', 'Fluxes', 1);
    obj.m.physics('NernstPlanckEquation').feature.create('BulkConcentrationsBC', 'Concentration', 1);
    obj.m.physics('NernstPlanckEquation').feature.create('continuity2dNernstPlanck', 'Continuity', 1);

    obj.m.physics.create('PoissonEquation', 'Electrostatics', 'geom');
    
    obj.m.physics('PoissonEquation').feature.create('chargeDensity', 'SpaceChargeDensity', 2);

    obj.m.physics('PoissonEquation').create('bpeChargeDensityBC', 'SurfaceChargeDensity', 1);        
    obj.m.physics('PoissonEquation').create('electrodePotentialBC', 'ElectricPotential', 1);  
    obj.m.physics('PoissonEquation').create('continuity2dPoisson', 'Continuity', 1);


    obj.m.physics.create('zeroSurfaceCurrent', 'GlobalEquations', 'geom');
    obj.m.physics('zeroSurfaceCurrent').create('zeroSurfaceCurrentEquation', 'GlobalEquations');
end