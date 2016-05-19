function update1dSelections(obj)
     %% define selections
    d = 1;
%     obj.surfaceNode.geom('geom',0); % point
%     obj.surfaceNode.set(1); % 1st created point as bpe surface
%     obj.surfaceNode.label('Surface node');
% 
%     obj.zetaNode.geom('geom',0); % point
%     obj.zetaNode.set(2); % 2nd created point as zeta plane
%     obj.zetaNode.label('Zeta plane node');
% 
%     obj.bulkNode.geom('geom',0); % point
%     obj.bulkNode.set(3); % 2rd created point far away at bulk solution
%     obj.bulkNode.label('Node at bulk');
    
    obj.allBoundaries.label('allBoundaries');
    obj.allBoundaries.geom('geom',d-1);
    obj.allBoundaries.all;

%     obj.bpeSurface.set('entitydim', d-1);
    obj.bpeSurface.label('bpeSurface');
    obj.bpeSurface.geom('geom',0); % point
    obj.bpeSurface.set(1); % 1st created point as bpe surface
    
    obj.entireSurface.set('entitydim', d-1);
    obj.entireSurface.set('input', {'bpeSurface'});
    obj.entireSurface.label('entireSurface');

    
%     obj.zetaPlane.set('entitydim', d-1);
    obj.zetaPlane.label('zetaPlane');
    obj.zetaPlane.geom('geom',0); % point
    obj.zetaPlane.set(2); % 2nd created point as zeta plane
    
    obj.bulkBoundary.set('entitydim', d-1);
    obj.bulkBoundary.label('bulkBoundary');
    obj.bulkBoundary.set('add', {'allBoundaries'});
    obj.bulkBoundary.set('subtract', {'entireSurface' 'zetaPlane'});
  
    %% interval selections
    obj.regionOfFirstDebyeLength.geom('geom',d); % interval
    obj.regionOfFirstDebyeLength.set(1); % ddl region of one debye length at surface
    obj.regionOfFirstDebyeLength.label('DDL region');
    
    obj.regionRemaining.geom('geom',d); 
    obj.regionRemaining.set(2);
    obj.regionRemaining.label('remaining region');
end