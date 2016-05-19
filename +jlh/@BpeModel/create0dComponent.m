function create0dComponent(obj)
    obj.m.modelNode.create('comp0d');
    obj.m.variable.create('var0d');
    obj.m.variable('var0d').model('comp0d');
    obj.m.physics.create('zeroNetCurrent0d', 'GlobalEquations');
end