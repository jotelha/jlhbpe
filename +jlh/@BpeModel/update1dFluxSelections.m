function update1dFluxSelections(obj)
     %% define selections
    d = 1;

    obj.allBoundaries.label('allBoundaries');
    obj.allBoundaries.geom('geom',d-1);
    obj.allBoundaries.all;

%     obj.entireSurface.set('entitydim', d-1);
%     obj.entireSurface.set('input', {'bpeSurface'});
%     obj.entireSurface.label('entireSurface');

    obj.workingElectrode.label('workingElectrode');
    obj.workingElectrode.geom('geom',d-1); % point
    obj.workingElectrode.set(1); 
    
    obj.counterElectrode.label('counterElectrode');
    obj.counterElectrode.geom('geom',d-1); % point
    obj.counterElectrode.set(2); 
    
    obj.lateralBoundary.label('lateralBoundary');
    obj.lateralBoundary.set('entitydim', d-1);
    obj.lateralBoundary.set('input', {'workingElectrode','counterElectrode'});
    
    obj.electrodes.label('electrodes');
    obj.electrodes.set('entitydim', d-1);
    obj.electrodes.set('input', {'workingElectrode','counterElectrode'});
end