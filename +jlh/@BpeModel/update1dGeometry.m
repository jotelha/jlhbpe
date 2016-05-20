function update1dGeometry(obj)
    fprintf('  Modifying geometry...\n');

    
%     obj.geom.repairTol(1e-16); % should be smaller than smallest elements
%     obj.space.set('intervals','many');
%     % space.set('p','0,lambdaD,L'); % an interval with three points
%     obj.space.set('p','0,epsilon,1'); % an interval with three points. dimensionless
% 
%     %electrode_surface = space.
% 
%     % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
%     obj.geom.run;
    obj.m.geom('geom1d').repairTol(1e-16); % should be smaller than smallest elements
%     obj.m.geom('geom1d').feature('space').set('intervals','many');
    % space.set('p','0,lambdaD,L'); % an interval with three points
%     obj.m.geom('geom1d').feature('space').set('p','0,epsilon,extendedDdlFactor*epsilon,1'); % an interval with three points. dimensionless
    obj.m.geom('geom1d').feature('ddl1d').set('p','0,epsilon'); % an interval with three points. dimensionless
    obj.m.geom('geom1d').feature('ddl1d').set('contributeto','ddl1dCumulative');
    
    obj.m.geom('geom1d').feature('extendedDdl1d').set('p','epsilon,extendedDdlFactor*epsilon'); % an interval with three points. dimensionless
    obj.m.geom('geom1d').feature('extendedDdl1d').set('contributeto','extendedDdl1dCumulative');
    
    obj.m.geom('geom1d').feature('space1d').set('p','extendedDdlFactor*epsilon,1'); % an interval with three points. dimensionless
    obj.m.geom('geom1d').feature('space1d').set('contributeto','space1dCumulative');
    
    % vertex selections
    obj.m.selection('surfaceVertex1d').label('surfaceVertex1d');
    obj.m.selection('surfaceVertex1d').geom('geom1d',0); % point
    obj.m.selection('surfaceVertex1d').set('entitydim', '0');
    obj.m.selection('surfaceVertex1d').set('xmin', '-epsilon/2');
    obj.m.selection('surfaceVertex1d').set('xmax', 'epsilon/2');
%     obj.m.selection('surfaceVertex1d').set(1); % 1st created point as bpe surface
    
    obj.m.selection('zetaVertex1d').label('zetaVertex1d');
    obj.m.selection('zetaVertex1d').geom('geom1d',0); % point
    obj.m.selection('zetaVertex1d').set('entitydim', '0');
    obj.m.selection('zetaVertex1d').set('xmin', 'epsilon/2');
    obj.m.selection('zetaVertex1d').set('xmax', '3*epsilon/2');
%     obj.m.selection('zetaVertex1d').set(2);
    
    obj.m.selection('bulkExitVertex1d').label('bulkExitVertex1d');
    obj.m.selection('bulkExitVertex1d').geom('geom1d',0); % point
    obj.m.selection('bulkExitVertex1d').set('entitydim', '0');
    obj.m.selection('bulkExitVertex1d').set('xmin', '1-epsilon/2');
    obj.m.selection('bulkExitVertex1d').set('xmax', '1+epsilon/2');
%     obj.m.selection('bulkExitVertex1d').set(3);
    % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
    
    % assembly 
    obj.m.geom('geom1d').feature('fin').set('action', 'assembly');
    obj.m.geom('geom1d').feature('fin').set('repairtol', '1.0E-16');
    obj.m.geom('geom1d').run('fin');
    
%     obj.m.geom('geom1d').run;
end