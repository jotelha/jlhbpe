function obj = create1dFluxPhysics(obj)
    d = 1;
     %% creating physics
    switch obj.transformation
        case 'logConcentration'
            fprintf('  Creating general PDE Nernst-Planck system for log concentration transform...\n');
            obj.NernstPlanckEquation = obj.m.physics.create('NernstPlanckEquation','GeneralFormPDE','geom');
        otherwise
            fprintf('  Creating coefficient PDE Nernst-Planck system...\n');
            obj.NernstPlanckEquation = obj.m.physics.create('NernstPlanckEquation','CoefficientFormPDE','geom');
    end

    % independent of transformation:
    obj.FluxAtElectrodes = obj.NernstPlanckEquation.create('FluxAtElectrodes', 'FluxBoundary',d-1);
    obj.BulkConcentrationsBC = obj.NernstPlanckEquation.create('BulkConcentrationsBC', 'DirichletBoundary', d-1);

    fprintf('  Creating coefficient PDE Poisson equation...\n');
    obj.PoissonEquation = obj.m.physics.create('PoissonEquation', 'CoefficientFormPDE', 'geom');
 
    obj.electrodePotentialBC = obj.PoissonEquation.create('electrodePotentialBC', 'DirichletBoundary', d-1);  

    obj.electrodeCurrent = obj.m.physics.create('electrodeCurrent', 'GlobalEquations', 'geom');
    obj.electrodeCurrentEquation = obj.electrodeCurrent.create('electrodeCurrentEquation', 'GlobalEquations');
end