function createFunctions(obj)
    %% create functions
    fprintf('  Creating global functions...\n');

    % this function can be used to manually impose initial values:
    obj.phiInterpolation = obj.m.func.create('phiInterpolation', 'Interpolation'); 
    obj.phiInterpolation.active(false); % deactivate per default

    % the analytical solution of Poisson-Boltzmann problem
    obj.phi_pb = obj.m.func.create('phi_pb', 'Analytic');
    obj.phi_pbx = obj.m.func.create('phi_pbx', 'Analytic');
    obj.c_pb = cell(obj.numberOfSpecies,1);
    obj.cInterpolation = cell(obj.numberOfSpecies,1); % aslo for initial value interpolation

    for i = 1:obj.numberOfSpecies   
        obj.c_pb{i} = obj.m.func.create(obj.c_pb_id{i}, 'Analytic');
        obj.cInterpolation{i} = obj.m.func.create(obj.cInterpolation_id{i}, 'Interpolation');
        obj.cInterpolation{i}.active(false); % deactivate per default
    end
    
    obj.m.func.create('smoothenBpeBC', 'Rectangle');
    
    obj.m.func.create('interpSolution1d', 'Interpolation');
    obj.m.func('interpSolution1d').active(false);
end