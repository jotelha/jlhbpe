function obj = create1dSelections(obj)
    fprintf('  Creating selections...\n');
    
    % point selections
%     obj.leftBoundaryOfSurface = obj.m.selection.create('leftBoundaryOfSurface', 'Box');
%     obj.rightBoundaryOfSurface = obj.m.selection.create('rightBoundaryOfSurface', 'Box');
%     obj.leftBoundaryOfZetaPlane = obj.m.selection.create('leftBoundaryOfZetaPlane', 'Box');
%     obj.rightBoundaryOfZetaPlane = obj.m.selection.create('rightBoundaryOfZetaPlane', 'Box');
%         
    obj.allBoundaries = obj.m.selection.create('allBoundaries', 'Explicit');
    obj.bpeSurface = obj.m.selection.create('bpeSurface', 'Explicit');

    obj.bulkBoundary = obj.m.selection.create('bulkBoundary', 'Difference');
%     obj.upperBoundary = obj.m.selection.create('upperBoundary', 'Difference');
    
%     obj.lateralBoundaryAtFirstDebyeLength = obj.m.selection.create('lateralBoundaryAtFirstDebyeLength', 'Box');
    
%     obj.lateralBoundaryAtRemainingDomain = obj.m.selection.create('lateralBoundaryAtRemainingDomain', 'Box');

%     obj.lateralBoundary = obj.m.selection.create('lateralBoundary', 'Union');

%     obj.openBoundary = obj.m.selection.create('openBoundary','Union');
    
    obj.zetaPlane = obj.m.selection.create('zetaPlane','Explicit'); % edge at some distance from BPE (one Debye length), representing the zeta plane

    
%     obj.electrodes = obj.m.selection.create('electrodes', 'Union');

    
%     obj.workingElectrode = obj.m.selection.create('workingElectrode', 'Box');
    
%     obj.counterElectrode = obj.m.selection.create('counterElectrode', 'Box');
    
%     obj.insulator = obj.m.selection.create('insulator', 'Adjacent');

%     obj.insulatorAdjacentToBpe = obj.m.selection.create('insulatorAdjacentToBpe', 'Adjacent');
        
    obj.entireSurface = obj.m.selection.create('entireSurface', 'Union');
    obj.reactingSurface = obj.m.selection.create('reactingSurface', 'Union');
   
    % domain selections
    obj.regionOfFirstDebyeLength = obj.m.selection.create('regionOfFirstDebyeLength','Explicit');
    obj.regionRemaining = obj.m.selection.create('regionRemaining','Explicit');

end