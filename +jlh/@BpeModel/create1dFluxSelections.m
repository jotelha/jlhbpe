function obj = create1dFluxSelections(obj)
    fprintf('  Creating selections...\n');
    
    % point selections
%     obj.leftBoundaryOfSurface = obj.m.selection.create('leftBoundaryOfSurface', 'Box');
%     obj.rightBoundaryOfSurface = obj.m.selection.create('rightBoundaryOfSurface', 'Box');
%     obj.leftBoundaryOfZetaPlane = obj.m.selection.create('leftBoundaryOfZetaPlane', 'Box');
%     obj.rightBoundaryOfZetaPlane = obj.m.selection.create('rightBoundaryOfZetaPlane', 'Box');
%         
    obj.allBoundaries = obj.m.selection.create('allBoundaries', 'Explicit');
    
    obj.workingElectrode = obj.m.selection.create('workingElectrode', 'Explicit');
    
    obj.counterElectrode = obj.m.selection.create('counterElectrode', 'Explicit');
    

    obj.lateralBoundary = obj.m.selection.create('lateralBoundary', 'Union');
    
    obj.electrodes = obj.m.selection.create('electrodes', 'Union'); 
end