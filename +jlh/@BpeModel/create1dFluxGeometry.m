function obj = create1dFluxGeometry(obj)
     %% creating geometry
    % geometry has to exist to define selections, etc
    fprintf('  Building geometry...\n');
    obj.geom = obj.m.geom.create('geom',1); % 1 for one dimension
    obj.space = obj.geom.create('space', 'Interval');
%     obj.ddl = obj.geom.create('ddl', 'Interval');

%     obj.leftEdgeOfBpe = obj.geom.create('leftEdgeOfBpe', 'Point');
%     obj.rightEdgeOfBpe = obj.geom.create('rightEdgeOfBpe', 'Point');
end