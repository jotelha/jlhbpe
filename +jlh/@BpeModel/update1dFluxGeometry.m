function obj = update1dFluxGeometry(obj)
    fprintf('  Modifying geometry...\n');

    obj.geom.repairTol(1e-16); % should be smaller than smallest elements
    obj.space.set('intervals','many');
    % space.set('p','0,lambdaD,L'); % an interval with three points
    obj.space.set('p','0,w'); % an interval with three points. dimensionless

    %electrode_surface = space.

    % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
    obj.geom.run;
end