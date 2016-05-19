function obj = create1dPhysics(obj)
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
    obj.FluxAtSurface = obj.NernstPlanckEquation.create('FluxAtSurface', 'FluxBoundary',d-1);
    obj.FluxAtElectrodes = obj.NernstPlanckEquation.create('FluxAtElectrodes', 'FluxBoundary',d-1);
    obj.FluxAtBulkBoundary = obj.NernstPlanckEquation.create('FluxAtBulkBoundary', 'FluxBoundary', d-1);
    obj.BulkConcentrationsBC = obj.NernstPlanckEquation.create('BulkConcentrationsBC', 'DirichletBoundary', d-1);

    fprintf('  Creating coefficient PDE Poisson equation...\n');
    obj.PoissonEquation = obj.m.physics.create('PoissonEquation', 'CoefficientFormPDE', 'geom');

    obj.bpeChargeDensityBC = obj.PoissonEquation.create('bpeChargeDensityBC', 'FluxBoundary', d-1);     
    obj.bpePotentialBC = obj.PoissonEquation.create('bpePotentialBC', 'DirichletBoundary', d-1);
    obj.entireSurfacePotentialGradientBC = obj.PoissonEquation.create('entireSurfacePotentialGradientBC', 'FluxBoundary', d-1);
    obj.insulatorSurfacePotentialGradientBC = obj.PoissonEquation.create('insulatorSurfacePotentialGradientBC', 'FluxBoundary', d-1);
    
    obj.bulkBoundaryPotentialGradientBC = obj.PoissonEquation.create('bulkBoundaryPotentialGradientBC', 'FluxBoundary', d-1);

    switch obj.bpePoissonEquationBC 
        case 'ChargeDensityBC'
            fprintf('    Robin boundary conditions are used...\n');
            % set surface charge at boundaries, Robin BC
            % phi + lambdaS * phi' = phi_bpe, or
            % C_Stern*(phi_bpe - phi) = eps0*eps* phi' = - sigma
            % with lambdaS width of Stern layer, sigma surface charge density
            obj.bpeChargeDensityBC.active(true);
            obj.bpePotentialBC.active(false);
        case 'PotentialBC'
            fprintf('    Dirichlet boundary conditions are used...\n');
            % potential at surface of bpe
            obj.bpeChargeDensityBC.active(false);
            obj.bpePotentialBC.active(true);
    end
    obj.bulkPotentialBC = obj.PoissonEquation.create('bulkPotentialBC', 'DirichletBoundary', d-1);
    obj.electrodePotentialBC = obj.PoissonEquation.create('electrodePotentialBC', 'DirichletBoundary', d-1);  

    obj.zeroSurfaceCurrent = obj.m.physics.create('zeroSurfaceCurrent', 'BoundaryODE', 'geom');
    obj.zeroSurfaceCurrentEquation = obj.zeroSurfaceCurrent.create('zeroSurfaceCurrentEquation', 'AlgebraicEquation', d-1);
end