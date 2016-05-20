function create1dComponent(obj)
    obj.m.modelNode.create('comp1d');
%     obj.m.variable.create('var0d');
%     obj.m.variable('var0d').model('comp0d');
%     obj.m.physics.create('zeroNetCurrent0d', 'GlobalEquations');

  %% creating geometry
    % geometry has to exist to define selections, etc
%     fprintf('  Building geometry...\n');
    obj.m.geom.create('geom1d',1); % 1 for one dimension
    obj.m.geom('geom1d').create('ddl1d', 'Interval');
    obj.m.geom('geom1d').create('extendedDdl1d', 'Interval');
    obj.m.geom('geom1d').create('space1d', 'Interval');
    
   % selections
    obj.m.geom('geom1d').selection.create('ddl1dCumulative', 'CumulativeSelection');
    obj.m.geom('geom1d').selection.create('extendedDdl1dCumulative', 'CumulativeSelection');
    obj.m.geom('geom1d').selection.create('space1dCumulative', 'CumulativeSelection');

    obj.m.selection.create('surfaceVertex1d', 'Box');
    obj.m.selection.create('zetaVertex1d', 'Box');
    obj.m.selection.create('bulkExitVertex1d', 'Box');

    % operators
    obj.m.cpl.create('projectReactionPlaneToSurface1d', 'LinearExtrusion', 'geom1d');
    obj.m.cpl.create('integrateSurface1d', 'Integration', 'geom1d');
    obj.m.cpl.create('integrateBulkExit1d', 'Integration', 'geom1d');

    % variables
    obj.m.variable.create('DomainVar1d');
    obj.m.variable('DomainVar1d').model('comp1d');

    obj.m.variable.create('SurfaceVar1d');
    obj.m.variable('SurfaceVar1d').model('comp1d');
    
    obj.m.variable.create('BulkVar1d');
    obj.m.variable('BulkVar1d').model('comp1d');

    % physics
    
    % weak formulation
    obj.m.physics.create('WeakFormulation1d','WeakFormPDE','geom1d');
    obj.m.physics('WeakFormulation1d').label('WeakFormulation1d');
    
    % continuities for assembly pairs
    obj.m.physics('WeakFormulation1d').feature.create('assemblyContinuity1d', 'Continuity', 0);


    % flux and concentration BC
    for i=1:obj.numberOfSpecies 
        obj.m.physics('WeakFormulation1d').create(obj.surfaceFluxBC_id{i}, 'WeakContribution', 0);
        obj.m.physics('WeakFormulation1d').feature(obj.surfaceFluxBC_id{i}).label(obj.surfaceFluxBC_id{i});
        
        obj.m.physics('WeakFormulation1d').create(obj.bulkFluxBC_id{i}, 'WeakContribution', 0);
        obj.m.physics('WeakFormulation1d').feature(obj.bulkFluxBC_id{i}).label(obj.bulkFluxBC_id{i});

        obj.m.physics('WeakFormulation1d').create(obj.bulkConcentrationBC_id{i}, 'WeakContribution', 0);
        obj.m.physics('WeakFormulation1d').feature(obj.bulkConcentrationBC_id{i}).label(obj.bulkConcentrationBC_id{i});
        obj.m.physics('WeakFormulation1d').feature(obj.bulkConcentrationBC_id{i}).create([obj.bulkConcentrationBC_id{i},'_aux'], 'AuxiliaryField', 0);
    end
    
    % Poisson BC
    obj.m.physics('WeakFormulation1d').create('BulkPotentialBC', 'WeakContribution', 0);
    obj.m.physics('WeakFormulation1d').feature('BulkPotentialBC').label('BulkPotentialBC');
    obj.m.physics('WeakFormulation1d').feature('BulkPotentialBC').create('BulkPotentialBC_aux', 'AuxiliaryField', 0);

%     obj.m.physics('WeakFormulation1d').create('ElectrodePotentialBC', 'WeakContribution', 1);
%     obj.m.physics('WeakFormulation1d').feature('ElectrodePotentialBC').label('ElectrodePotentialBC');
%     obj.m.physics('WeakFormulation1d').feature('ElectrodePotentialBC').create('ElectrodePotentialBC_aux', 'AuxiliaryField', 1);
 
    obj.m.physics('WeakFormulation1d').create('SurfaceChargeDensityBC', 'WeakContribution', 0);
    obj.m.physics('WeakFormulation1d').feature('SurfaceChargeDensityBC').label('SurfaceChargeDensityBC');
    
    
    % classical formulation
%     obj.m.physics.create('NernstPlanckEquation1d','CoefficientFormPDE','geom1d');
%     obj.m.physics('NernstPlanckEquation1d').create('FluxAtSurface1d', 'FluxBoundary',0);
%     obj.m.physics('NernstPlanckEquation1d').create('BulkConcentrationsBC1d', 'DirichletBoundary', 0);
% 
%     obj.m.physics.create('PoissonEquation1d', 'CoefficientFormPDE', 'geom1d');
%     obj.m.physics('PoissonEquation1d').create('bpeChargeDensityBC1d', 'FluxBoundary', 0);     
%     obj.m.physics('PoissonEquation1d').create('bpePotentialBC1d', 'DirichletBoundary', 0);
%     obj.m.physics('PoissonEquation1d').feature('bpePotentialBC1d').active(false);
%     obj.m.physics('PoissonEquation1d').create('bulkPotentialBC1d', 'DirichletBoundary', 0);
   
    obj.m.physics.create('zeroNetCurrent1d', 'BoundaryODE', 'geom1d');
    obj.m.physics('zeroNetCurrent1d').create('zeroSurfaceCurrentEquation1d', 'AlgebraicEquation', 0);
   
    obj.m.mesh.create('standardMesh1d', 'geom1d');    
    obj.m.mesh('standardMesh1d').create('meshingEdge1d', 'Edge'); 
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').create('meshDistributionDDL1d', 'Distribution');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').create('meshDistributionExtendedDDL1d', 'Distribution');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').create('meshDistributionRemaining1d', 'Distribution');
    
    obj.m.mesh.create('explicitMesh1d', 'geom1d');    
    obj.m.mesh('explicitMesh1d').create('meshingEdge1d', 'Edge'); 
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').create('meshDistributionDDL1d', 'Distribution');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').create('meshDistributionExtendedDDL1d', 'Distribution');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').create('meshDistributionRemaining1d', 'Distribution');

end