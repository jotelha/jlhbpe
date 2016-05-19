function obj = createGeometry(obj)
     %% creating geometry
    % geometry has to exist to define selections, etc
    fprintf('  Building geometry...\n');
    obj.geom = obj.m.geom.create('geom',2); % 1 for one dimension
    % space = geom.feature.create('space','Interval');
    % obj.space = obj.geom.feature.create('space', 'Rectangle');
    obj.space = obj.geom.create('space', 'Rectangle');
    obj.ddl = obj.geom.create('ddl', 'Rectangle');
    
    obj.gapLeft = obj.geom.create('gapLeft', 'Rectangle');
    obj.gapRight = obj.geom.create('gapRight', 'Rectangle');
    
    obj.spaceArray = obj.geom.create('spaceArray','Array');


    obj.leftEdgeOfBpe = obj.geom.create('leftEdgeOfBpe', 'Point');
    obj.rightEdgeOfBpe = obj.geom.create('rightEdgeOfBpe', 'Point');
    
    %% for explicit electrodes
    obj.leftEdgeOfCathode = obj.geom.create('leftEdgeOfCathode', 'Point');
    obj.rightEdgeOfCathode = obj.geom.create('rightEdgeOfCathode', 'Point');
    
    obj.leftEdgeOfAnode = obj.geom.create('leftEdgeOfAnode', 'Point');
    obj.rightEdgeOfAnode = obj.geom.create('rightEdgeOfAnode', 'Point');

    % space.set('intervals','many');
    % space.set('p','0,lambdaD,L'); % an interval with three points
    % % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
    % geom.run;
end