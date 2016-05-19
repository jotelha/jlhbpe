function obj = updateDatasetsForChoppedGeometry(obj,dset)
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
    
    obj.bpeSurfaceResults.set('data', dset);
    obj.bpeSurfaceResults.label('bpeSurfaceResults');
    obj.bpeSurfaceResults.selection.named('geom_bpeSurface');
   
    obj.entireSurfaceResults.set('data', dset);
    obj.entireSurfaceResults.label('entireSurfaceResults');
	obj.entireSurfaceResults.selection.named('geom_entireSurface');

    obj.zetaPlaneResults.set('data', dset);
    obj.zetaPlaneResults.label('zetaPlaneResults');
	obj.zetaPlaneResults.selection.named('geom_zetaPlane');
    
    obj.bulkBoundaryResults.set('data', dset);
    obj.bulkBoundaryResults.label('bulkBoundaryResults');
% 	obj.bulkBoundaryResults.selection.named('bulkBoundary');
% fix for plots 2016-04-07
	obj.bulkBoundaryResults.selection.named('geom_upperBoundary');
    
    obj.weResults.set('data', dset);
    obj.weResults.label('weResults');
	obj.weResults.selection.named('geom_workingElectrode');

    obj.ceResults.set('data', dset);
    obj.ceResults.label('ceResults');
    obj.ceResults.selection.named('geom_counterElectrode');

    obj.centralCrossectionResults.set('data', dset);
    obj.centralCrossectionResults.label('centralCrossectionResults');
    obj.centralCrossectionResults.set('genpoints', {'0' '0'; '0' '1'});
    obj.centralCrossectionResults.set('spacevars', 'cx');
 
    obj.centralDDLCrossectionResults.set('data', dset);
    obj.centralDDLCrossectionResults.label('centralDDLCrossectionResults');
    obj.centralDDLCrossectionResults.set('genpoints', {'0' '0'; '0' num2str(obj.epsilon)});
    obj.centralDDLCrossectionResults.set('spacevars', 'cx');
    
    obj.leftBpeEdgeCrossectionResults.set('data', dset);
    obj.leftBpeEdgeCrossectionResults.label('leftBpeEdgeCrossectionResults');
    obj.leftBpeEdgeCrossectionResults .set('genpoints', {'-w_bpe/2' '0'; '-w_bpe/2' '1'});
    obj.leftBpeEdgeCrossectionResults.set('spacevars', 'cx');
    
    obj.rightBpeEdgeCrossectionResults.set('data', dset);
    obj.rightBpeEdgeCrossectionResults.label('rightBpeEdgeCrossectionResults');
    obj.rightBpeEdgeCrossectionResults .set('genpoints', {'w_bpe/2' '0'; 'w_bpe/2' '1'});
    obj.rightBpeEdgeCrossectionResults.set('spacevars', 'cx');
    
    obj.cathodeCrossectionResults.set('data', dset);
    obj.cathodeCrossectionResults.label('cathodeCrossectionResults');
    obj.cathodeCrossectionResults .set('genpoints', {'-w_bpe/2-w_insulatorLeft-w_cathode/2' '0'; '-w_bpe/2-w_insulatorLeft-w_cathode/2' '1'});
    obj.cathodeCrossectionResults.set('spacevars', 'cx');obj.cathodeCrossectionResults.set('data', dset);
    
    obj.anodeCrossectionResults.set('data', dset);
    obj.anodeCrossectionResults.label('anodeCrossectionResults');
    obj.anodeCrossectionResults .set('genpoints', {'w_bpe/2+w_insulatorRight+w_anode/2' '0'; 'w_bpe/2+w_insulatorRight+w_anode/2' '1'});
    obj.anodeCrossectionResults.set('spacevars', 'cx');
    
    obj.halfDebyeLengthResults.set('data', dset);
    obj.halfDebyeLengthResults.label('halfDebyeLengthResults');
    obj.halfDebyeLengthResults .set('genpoints', {'-w/2' 'epsilon/2'; 'w/2' 'epsilon/2'});
    obj.halfDebyeLengthResults.set('spacevars', 'cx');
end