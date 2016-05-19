function obj = createSelections(obj)
    fprintf('  Creating selections...\n');
%     obj.surfaceNode = obj.m.selection.create('surfaceNode','Explicit'); % point at surface of bpe
%     obj.bulkNode = obj.m.selection.create('bulkNode','Explicit'); % point at bulk

    % in order to allow plot preparation (especially to enable duplicating
    % plots with assigned selections):
    % surfaceNode.geom('geom',0); % point
    % surfaceNode.set(1); % 1st created point as bpe surface
    % surfaceNode.label('Surface node');
    % 
    % zetaNode.geom('geom',0); % point
    % zetaNode.set(2); % 2nd created point as zeta plane
    % zetaNode.label('Zeta plane node');
    % 
    % bulkNode.geom('geom',0); % point
    % bulkNode.set(3); % 2rd created point far away at bulk solution
    % bulkNode.label('Node at bulk');
    % 
    
    % point selections
    obj.leftBoundaryOfSurface = obj.m.selection.create('leftBoundaryOfSurface', 'Box');
    obj.rightBoundaryOfSurface = obj.m.selection.create('rightBoundaryOfSurface', 'Box');
    obj.leftBoundaryOfZetaPlane = obj.m.selection.create('leftBoundaryOfZetaPlane', 'Box');
    obj.rightBoundaryOfZetaPlane = obj.m.selection.create('rightBoundaryOfZetaPlane', 'Box');
    
    % edge selections
    
    obj.allBoundaries = obj.m.selection.create('allBoundaries', 'Explicit');
    
    obj.bpeSurface = obj.m.selection.create('bpeSurface', 'Box');
%     obj.bpeSurface.set('entitydim', '1');
%     obj.bpeSurface.label('bpeSurface');
%     obj.bpeSurface.set('xmin', '0');
%     obj.bpeSurface.set('xmax', '0');
%     obj.bpeSurface.set('ymax', '0');
%     obj.bpeSurface.set('ymax', '0');
    
%     obj.bulkBoundary = obj.m.selection.create('bulkBoundary', 'Box');
    obj.bulkBoundary = obj.m.selection.create('bulkBoundary', 'Difference');
    obj.upperBoundary = obj.m.selection.create('upperBoundary', 'Difference');
    
    obj.lateralBoundaryAtFirstDebyeLength = obj.m.selection.create('lateralBoundaryAtFirstDebyeLength', 'Box');
%     obj.electrodesAtFirstDebyeLength.set('entitydim', '1');
%     obj.electrodesAtFirstDebyeLength.label('electrodes at first debye length');
%     obj.electrodesAtFirstDebyeLength.set('ymin', 'epsilon/2');
%     obj.electrodesAtFirstDebyeLength.set('ymax', 'epsilon/2');
    
    obj.lateralBoundaryAtRemainingDomain = obj.m.selection.create('lateralBoundaryAtRemainingDomain', 'Box');
%     obj.electrodesAtRemainingDomain.set('entitydim', '1');
%     obj.electrodesAtRemainingDomain.label('electrodes at remianing domain');
%     obj.electrodesAtRemainingDomain.set('ymin', '1/2');
%     obj.electrodesAtRemainingDomain.set('ymax', '1/2');
%     obj.bulkBoundary.set('entitydim', '1');
%     obj.bulkBoundary.label('bulkBoundary');
%     obj.bulkBoundary.set('xmin', '0');
%     obj.bulkBoundary.set('xmax', '0');
%     obj.bulkBoundary.set('ymin', '1');
%     obj.bulkBoundary.set('ymax', '1');
    obj.lateralBoundary = obj.m.selection.create('lateralBoundary', 'Union');

%     obj.openBoundary = obj.m.selection.create('openBoundary','Union');
    
    obj.zetaPlane = obj.m.selection.create('zetaPlane','Box'); % edge at some distance from BPE (one Debye length), representing the zeta plane
%     obj.zetaPlane.set('entitydim', '1');
%     obj.zetaPlane.label('zetaPlane');
%     obj.zetaPlane.set('xmin', '0');
%     obj.zetaPlane.set('xmax', '0');
%     obj.zetaPlane.set('ymin', 'epsilon');
%     obj.zetaPlane.set('ymax', 'epsilon');
    
 
    
    obj.electrodes = obj.m.selection.create('electrodes', 'Union');
%     obj.electrodes.set('entitydim', '1');
%     obj.electrodes.set('input', {'electrodesAtFirstDebyeLength' 'electrodesAtRemainingDomain'});
%     obj.electrodes.label('electrodes');
    
    obj.workingElectrode = obj.m.selection.create('workingElectrode', 'Box');
%     obj.workingElectrode.set('entitydim', '1');
%     obj.workingElectrode.label('workingElectrode');
%     obj.workingElectrode.set('ymin', '1/2');
%     obj.workingElectrode.set('ymax', '1/2');
%     obj.workingElectrode.set('xmax', '0');
    
    obj.counterElectrode = obj.m.selection.create('counterElectrode', 'Box');
%     obj.counterElectrode.set('entitydim', '1');
%     obj.counterElectrode.label('counterElectrode');
%     obj.counterElectrode.set('ymin', '1/2');
%     obj.counterElectrode.set('ymax', '1/2');
%     obj.counterElectrode.set('xmin', '0');
    
    obj.insulator = obj.m.selection.create('insulator', 'Adjacent');

    obj.insulatorAdjacentToBpe = obj.m.selection.create('insulatorAdjacentToBpe', 'Adjacent');
%     obj.insulatorAdjacentToBpe.set('entitydim', '1');
%     obj.insulatorAdjacentToBpe.set('input', {'bpeSurface'});
%     obj.insulatorAdjacentToBpe.label('insulatorAdjacentToBpe');
        
    obj.entireSurface = obj.m.selection.create('entireSurface', 'Union');
%     obj.entireSurface.set('entitydim', '1');
%     obj.entireSurface.set('input', {'insulatorAdjacentToBpe' 'bpeSurface'});
%     obj.entireSurface.label('entireSurface');
    
   obj.reactingSurface = obj.m.selection.create('reactingSurface', 'Union');
   
    % domain selections
    obj.regionOfFirstDebyeLength = obj.m.selection.create('regionOfFirstDebyeLength','Adjacent');
%     obj.regionOfFirstDebyeLength.set('entitydim', '2');
%     % obj.regionOfFirstDebyeLength.set('input', {'bpeSurface'});
% %     obj.regionOfFirstDebyeLength.geom('geom',2); % domain
% %     obj.regionOfFirstDebyeLength.set(1); % ddl region of one debye length at surface
%     obj.regionOfFirstDebyeLength.label('DDL region');
end