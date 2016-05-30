function obj = createPhysicsForWeakForm(obj)
    obj.m.physics.create('WeakFormulation','WeakFormPDE','geom');
    obj.m.physics('WeakFormulation').label('WeakFormulation');

    
    % flux and concentration BC
    for i=1:obj.numberOfSpecies 
        obj.m.physics('WeakFormulation').create(obj.surfaceFluxBC_id{i}, 'WeakContribution', 1);
        obj.m.physics('WeakFormulation').feature(obj.surfaceFluxBC_id{i}).label(obj.surfaceFluxBC_id{i});

        obj.m.physics('WeakFormulation').create(obj.bulkConcentrationBC_id{i}, 'PointwiseConstraint', 1);
        obj.m.physics('WeakFormulation').feature(obj.bulkConcentrationBC_id{i}).label(obj.bulkConcentrationBC_id{i});
%         obj.m.physics('WeakFormulation').feature(obj.bulkConcentrationBC_id{i}).create([obj.bulkConcentrationBC_id{i},'_aux'], 'AuxiliaryField', 1);
    end
    
    % Poisson BC
%     obj.m.physics('WeakFormulation').create('BulkPotentialBC', 'WeakContribution', 1);
%     obj.m.physics('WeakFormulation').feature('BulkPotentialBC').label('BulkPotentialBC');
%     obj.m.physics('WeakFormulation').feature('BulkPotentialBC').create('BulkPotentialBC_aux', 'AuxiliaryField', 1);

    obj.m.physics('WeakFormulation').create('ElectrodePotentialBC', 'PointwiseConstraint', 1);
    obj.m.physics('WeakFormulation').feature('ElectrodePotentialBC').label('ElectrodePotentialBC');
%     obj.m.physics('WeakFormulation').feature('ElectrodePotentialBC').create('ElectrodePotentialBC_aux', 'AuxiliaryField', 1);
    
    
    obj.m.physics('WeakFormulation').create('SurfaceChargeDensityBC', 'WeakContribution', 1);
    obj.m.physics('WeakFormulation').feature('SurfaceChargeDensityBC').label('SurfaceChargeDensityBC');

    % Identity pair continuity
    obj.m.physics('WeakFormulation').create('zetaPlaneContinuity', 'Continuity', 1);

    obj.zeroSurfaceCurrent = obj.m.physics.create('zeroSurfaceCurrent', 'GlobalEquations', 'geom');
    obj.zeroSurfaceCurrentEquation = obj.zeroSurfaceCurrent.create('zeroSurfaceCurrentEquation', 'GlobalEquations');
end