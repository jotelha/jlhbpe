function obj = updateSelections(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    %% define selections
%     obj.surfaceNode.geom('geom',0); % point
%     obj.surfaceNode.set(1); % 1st created point as bpe surface
%     obj.surfaceNode.label('Surface node');
% 
%     obj.zetaNode.geom('geom',0); % point
%     obj.zetaNode.set(2); % 2nd created point as zeta plane
%     obj.zetaNode.label('Zeta plane node');
% 
%     obj.bulkBoundary.geom('geom',0); % point
%     obj.bulkBoundary.set(3); % 2rd created point far away at bulk solution
%     obj.bulkBoundary.label('Node at bulk');
% 
%     obj.regionOfFirstDebyeLength.geom('geom',1); % interval
%     obj.regionOfFirstDebyeLength.set(1); % ddl region of one debye length at surface
%     obj.regionOfFirstDebyeLength.label('DDL region');

    % point selections
        
    obj.leftBoundaryOfSurface.set('entitydim', '0');
    obj.leftBoundaryOfSurface.label('leftBoundaryOfSurface');
    obj.leftBoundaryOfSurface.set('xmin', '-w/2');
    obj.leftBoundaryOfSurface.set('xmax', '-w/2');
    obj.leftBoundaryOfSurface.set('ymin', '0');
    obj.leftBoundaryOfSurface.set('ymax', '0');
    
    obj.rightBoundaryOfSurface.set('entitydim', '0');
    obj.rightBoundaryOfSurface.label('rightBoundaryOfSurface');
    obj.rightBoundaryOfSurface.set('xmin', 'w/2');
    obj.rightBoundaryOfSurface.set('xmax', 'w/2');
    obj.rightBoundaryOfSurface.set('ymin', '0');
    obj.rightBoundaryOfSurface.set('ymax', '0');
    
    obj.leftBoundaryOfZetaPlane.set('entitydim', '0');
    obj.leftBoundaryOfZetaPlane.label('leftBoundaryOfZetaPlane');
    obj.leftBoundaryOfZetaPlane.set('xmin', '-w/2');
    obj.leftBoundaryOfZetaPlane.set('xmax', '-w/2');
    obj.leftBoundaryOfZetaPlane.set('ymin', 'epsilon');
    obj.leftBoundaryOfZetaPlane.set('ymax', 'epsilon');
    
    obj.rightBoundaryOfZetaPlane.set('entitydim', '0');
    obj.rightBoundaryOfZetaPlane.label('rightBoundaryOfZetaPlane');
    obj.rightBoundaryOfZetaPlane.set('xmin', 'w/2');
    obj.rightBoundaryOfZetaPlane.set('xmax', 'w/2');
    obj.rightBoundaryOfZetaPlane.set('ymin', 'epsilon');
    obj.rightBoundaryOfZetaPlane.set('ymax', 'epsilon');
    
%     % edge selections
%     obj.allBoundaries.label('allBoundaries');
%     obj.allBoundaries.geom('geom',1);
%     obj.allBoundaries.all;
% 
%     obj.bpeSurface.set('entitydim', '1');
%     obj.bpeSurface.label('bpeSurface');
%     obj.bpeSurface.set('xmin', '0');
%     obj.bpeSurface.set('xmax', '0');
%     obj.bpeSurface.set('ymin', '0');
%     obj.bpeSurface.set('ymax', '0');
%     
%     obj.insulatorAdjacentToBpe.set('entitydim', '1');
%     obj.insulatorAdjacentToBpe.set('input', {'bpeSurface'});
%     obj.insulatorAdjacentToBpe.label('insulatorAdjacentToBpe');
%    
%     %% zeta plane
%     
%     obj.zetaPlane.set('entitydim', '1');
%     obj.zetaPlane.label('zetaPlane');
%     obj.zetaPlane.set('xmin', '0');
%     obj.zetaPlane.set('xmax', '0');
%     obj.zetaPlane.set('ymin', '3/4*epsilon');
%     obj.zetaPlane.set('ymax', '5/4*epsilon');
%     % strangely, factors are necessary, otherwise selection somtetimes
%     % empty
%     
%     %% lateral boundaries
%       obj.lateralBoundaryAtFirstDebyeLength.set('entitydim', '1');
%     obj.lateralBoundaryAtFirstDebyeLength.label('lateral boundary at first debye length');
%     obj.lateralBoundaryAtFirstDebyeLength.set('ymin', 'epsilon/2');
%     obj.lateralBoundaryAtFirstDebyeLength.set('ymax', 'epsilon/2');
% 
%     obj.lateralBoundaryAtRemainingDomain.set('entitydim', '1');
%     obj.lateralBoundaryAtRemainingDomain.label('lateral boundary at remaining domain');
%     obj.lateralBoundaryAtRemainingDomain.set('ymin', '1/2');
%     obj.lateralBoundaryAtRemainingDomain.set('ymax', '1/2');
%     
%     obj.lateralBoundary.set('entitydim', '1');
%     obj.lateralBoundary.set('input', {'lateralBoundaryAtFirstDebyeLength' 'lateralBoundaryAtRemainingDomain'});
%     obj.lateralBoundary.label('lateralBoundary');
%    
%     if(~obj.explicitElectrodeGeometry)
%         obj.reactingSurface.set('entitydim', '1');
%         obj.reactingSurface.set('input', {'bpeSurface'});
%         obj.reactingSurface.label('reactingSurface');
%         
% %         obj.bulkBoundary.set('entitydim', '1');
% %         obj.bulkBoundary.label('bulkBoundary');
% %         obj.bulkBoundary.set('xmin', '0');
% %         obj.bulkBoundary.set('xmax', '0');
% %         obj.bulkBoundary.set('ymin', '1');
% %         obj.bulkBoundary.set('ymax', '1');
%     
% 
%         obj.electrodes.set('entitydim', '1');
%         obj.electrodes.set('input', {'lateralBoundaryAtFirstDebyeLength' 'lateralBoundaryAtRemainingDomain'});
%         obj.electrodes.label('electrodes');
% 
%         obj.workingElectrode.set('entitydim', '1');
%         obj.workingElectrode.label('workingElectrode');
%         obj.workingElectrode.set('condition', 'inside');
%         obj.workingElectrode.set('xmin', '-Inf');
%         obj.workingElectrode.set('xmax', '-w/2');
%         obj.workingElectrode.set('ymin', '-Inf');
%         obj.workingElectrode.set('ymax', 'Inf');
% 
%         obj.counterElectrode.set('entitydim', '1');
%         obj.counterElectrode.label('counterElectrode');
%         obj.counterElectrode.set('condition', 'inside');
%         obj.counterElectrode.set('xmin', 'w/2');
%         obj.counterElectrode.set('xmax', 'Inf');
%         obj.counterElectrode.set('ymin', '-Inf');
%         obj.counterElectrode.set('ymax', 'Inf');
%         
%         obj.insulator.set('entitydim', '1');
%         obj.insulator.set('input', {'bpeSurface'});
%         obj.insulator.label('insulator');
%         
%         obj.entireSurface.set('entitydim', '1');
%         obj.entireSurface.set('input', {'insulator' 'bpeSurface'});
%         obj.entireSurface.label('entireSurface');
%     else   
%         obj.workingElectrode.set('entitydim', '1');
%         obj.workingElectrode.label('workingElectrode');
% %         obj.workingElectrode.set('condition', 'inside');
%         obj.workingElectrode.set('xmin', '-w_bpe/2-w_insulatorLeft-w_cathode/2');
%         obj.workingElectrode.set('xmax', '-w_bpe/2-w_insulatorLeft-w_cathode/2');
%         obj.workingElectrode.set('ymin', '0');
%         obj.workingElectrode.set('ymax', '0');
% 
%         obj.counterElectrode.set('entitydim', '1');
%         obj.counterElectrode.label('counterElectrode');
% %         obj.counterElectrode.set('condition', 'inside');
%         obj.counterElectrode.set('xmin', 'w_bpe/2+w_insulatorRight+w_anode/2');
%         obj.counterElectrode.set('xmax', 'w_bpe/2+w_insulatorRight+w_anode/2');
%         obj.counterElectrode.set('ymin', '0');
%         obj.counterElectrode.set('ymax', '0');
%         
%         obj.electrodes.set('entitydim', '1');
%         obj.electrodes.set('input', {'counterElectrode' 'workingElectrode'});
%         obj.electrodes.label('electrodes');
%         
%         obj.reactingSurface.set('entitydim', '1');
%         obj.reactingSurface.set('input', {'bpeSurface' 'counterElectrode' 'workingElectrode'});
%         obj.reactingSurface.label('reactingSurface');
%         
%         obj.insulator.set('entitydim', '1');
%         obj.insulator.set('input', {'bpeSurface' 'electrodes'});
%         obj.insulator.label('insulator');
%         
%         obj.entireSurface.set('entitydim', '1');
%         obj.entireSurface.set('input', {'insulator' 'bpeSurface' 'electrodes'});
%         obj.entireSurface.label('entireSurface');
%     end
%     
%     obj.bulkBoundary.set('entitydim', '1');
%     obj.bulkBoundary.label('bulkBoundary');
% %         obj.bulkBoundary.set('xmin', '-Inf');
% %         obj.bulkBoundary.set('xmax', 'Inf');
% %         obj.bulkBoundary.set('ymin', 'epsilon/2');
% %         obj.bulkBoundary.set('ymax', 'Inf');
%     obj.bulkBoundary.set('add', {'allBoundaries'});
%     obj.bulkBoundary.set('subtract', {'entireSurface' 'zetaPlane'});
%     
%     obj.upperBoundary.set('entitydim', '1');
%     obj.upperBoundary.label('upperBoundary');
%     obj.upperBoundary.set('add', {'bulkBoundary'});
%     obj.upperBoundary.set('subtract', {'lateralBoundary'});
%     
%     % domain selections
%     obj.regionOfFirstDebyeLength.set('entitydim', '1');
%     obj.regionOfFirstDebyeLength.set('outputdim', '2');
%     obj.regionOfFirstDebyeLength.set('input', {'bpeSurface'});
%     obj.regionOfFirstDebyeLength.label('DDL region');
end