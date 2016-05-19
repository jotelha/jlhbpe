function createVariables(obj)
    obj.domainVariables = obj.m.variable.create('domainVariables');
    obj.surfaceVariables = obj.m.variable.create('surfaceVariables');
    obj.bulkVariables = obj.m.variable.create('bulkVariables');
    obj.upperBoundaryVariables = obj.m.variable.create('upperBoundaryVariables');


    obj.electrodeVariables = obj.m.variable.create('electrodeVariables');
    obj.weVariables = obj.m.variable.create('weVariables');
    obj.ceVariables = obj.m.variable.create('ceVariables');
end