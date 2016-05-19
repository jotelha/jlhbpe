function obj = updateDatasets(obj,dset)
%     obj.weResults                   = obj.m.result.dataset.create('weResults', 'Edge2D');
%     obj.ceResults                   = obj.m.result.dataset.create('ceResults', 'Edge2D');
%     obj.bpeSurfaceResults           = obj.m.result.dataset.create('bpeSurfaceResults', 'Edge2D');
%     obj.zetaPlaneResults            = obj.m.result.dataset.create('zetaPlaneResults', 'Edge2D');
%     obj.bulkBoundaryResults         = obj.m.result.dataset.create('bulkBoundaryResults', 'Edge2D');    
%     obj.entireSurfaceResults        = obj.m.result.dataset.create('entireSurfaceResults', 'Edge2D');
% 	obj.centralCrossectionResults   = obj.m.result.dataset.create('centralCrossectionResults', 'CutLine2D');
%     

    if( nargin < 2 )
        dset = obj.getLatestDataset();
    end
    
    obj.m.result.dataset('bpeSurfaceResults').set('data', dset);
    obj.m.result.dataset('bpeSurfaceResults').label('bpeSurfaceResults');
    obj.m.result.dataset('bpeSurfaceResults').selection.named('geom_bpeSurface');
   
    obj.m.result.dataset('entireSurfaceResults').set('data', dset);
    obj.m.result.dataset('entireSurfaceResults').label('entireSurfaceResults');
	obj.m.result.dataset('entireSurfaceResults').selection.named('geom_entireSurface');

    obj.m.result.dataset('zetaPlaneResults').set('data', dset);
    obj.m.result.dataset('zetaPlaneResults').label('zetaPlaneResults');
	obj.m.result.dataset('zetaPlaneResults').selection.named('geom_zetaPlane');
    
    obj.m.result.dataset('bulkBoundaryResults').set('data', dset);
    obj.m.result.dataset('bulkBoundaryResults').label('bulkBoundaryResults');
% 	obj.m.result.dataset('bulkBoundaryResults.selection.named('bulkBoundary');
% fix for plots 2016-04-07
	obj.m.result.dataset('bulkBoundaryResults').selection.named('geom_upperBoundary');
    
    obj.m.result.dataset('weResults').set('data', dset);
    obj.m.result.dataset('weResults').label('weResults');
	obj.m.result.dataset('weResults').selection.named('geom_workingElectrode');

    obj.m.result.dataset('ceResults').set('data', dset);
    obj.m.result.dataset('ceResults').label('ceResults');
    obj.m.result.dataset('ceResults').selection.named('geom_counterElectrode');

    obj.m.result.dataset('centralCrossectionResults').set('data', dset);
    obj.m.result.dataset('centralCrossectionResults').label('centralCrossectionResults');
    obj.m.result.dataset('centralCrossectionResults').set('genpoints', {'0' '0'; '0' '1'});
    obj.m.result.dataset('centralCrossectionResults').set('spacevars', 'cx');
 
    obj.m.result.dataset('centralDDLCrossectionResults').set('data', dset);
    obj.m.result.dataset('centralDDLCrossectionResults').label('centralDDLCrossectionResults');
    obj.m.result.dataset('centralDDLCrossectionResults').set('genpoints', {'0' '0'; '0' num2str(obj.epsilon)});
    obj.m.result.dataset('centralDDLCrossectionResults').set('spacevars', 'cx');
    
    obj.m.result.dataset('leftBpeEdgeCrossectionResults').set('data', dset);
    obj.m.result.dataset('leftBpeEdgeCrossectionResults').label('leftBpeEdgeCrossectionResults');
    obj.m.result.dataset('leftBpeEdgeCrossectionResults').set('genpoints', {'-w_bpe/2' '0'; '-w_bpe/2' '1'});
    obj.m.result.dataset('leftBpeEdgeCrossectionResults').set('spacevars', 'cx');
    
    obj.m.result.dataset('rightBpeEdgeCrossectionResults').set('data', dset);
    obj.m.result.dataset('rightBpeEdgeCrossectionResults').label('rightBpeEdgeCrossectionResults');
    obj.m.result.dataset('rightBpeEdgeCrossectionResults').set('genpoints', {'w_bpe/2' '0'; 'w_bpe/2' '1'});
    obj.m.result.dataset('rightBpeEdgeCrossectionResults').set('spacevars', 'cx');
    
    obj.m.result.dataset('cathodeCrossectionResults').set('data', dset);
    obj.m.result.dataset('cathodeCrossectionResults').label('cathodeCrossectionResults');
    obj.m.result.dataset('cathodeCrossectionResults').set('genpoints', {'-w_bpe/2-w_insulatorLeft-w_cathode/2' '0'; '-w_bpe/2-w_insulatorLeft-w_cathode/2' '1'});
    obj.m.result.dataset('cathodeCrossectionResults').set('spacevars', 'cx');
%     obj.m.result.dataset('cathodeCrossectionResults').set('data', dset);
    
    obj.m.result.dataset('anodeCrossectionResults').set('data', dset);
    obj.m.result.dataset('anodeCrossectionResults').label('anodeCrossectionResults');
    obj.m.result.dataset('anodeCrossectionResults').set('genpoints', {'w_bpe/2+w_insulatorRight+w_anode/2' '0'; 'w_bpe/2+w_insulatorRight+w_anode/2' '1'});
    obj.m.result.dataset('anodeCrossectionResults').set('spacevars', 'cx');
    
    obj.m.result.dataset('halfDebyeLengthResults').set('data', dset);
    obj.m.result.dataset('halfDebyeLengthResults').label('halfDebyeLengthResults');
    obj.m.result.dataset('halfDebyeLengthResults').set('genpoints', {'-w/2' 'epsilon/2'; 'w/2' 'epsilon/2'});
    obj.m.result.dataset('halfDebyeLengthResults').set('spacevars', 'cx');
end