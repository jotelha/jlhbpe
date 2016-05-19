function obj = prepareDatasets(obj)
%     obj.weResults                   = obj.m.result.dataset.create('weResults', 'Edge2D');
%     obj.ceResults                   = obj.m.result.dataset.create('ceResults', 'Edge2D');
%     obj.bpeSurfaceResults           = obj.m.result.dataset.create('bpeSurfaceResults', 'Edge2D');
%     obj.zetaPlaneResults            = obj.m.result.dataset.create('zetaPlaneResults', 'Edge2D');
%     obj.bulkBoundaryResults         = obj.m.result.dataset.create('bulkBoundaryResults', 'Edge2D');    
%     obj.entireSurfaceResults        = obj.m.result.dataset.create('entireSurfaceResults', 'Edge2D');
% 	obj.centralCrossectionResults   = obj.m.result.dataset.create('centralCrossectionResults', 'CutLine2D');
%     
    obj.bpeSurfaceResults.label('bpeSurfaceResults');
    obj.bpeSurfaceResults.selection.named('bpeSurface');
   
    obj.entireSurfaceResults.label('entireSurfaceResults');
	obj.entireSurfaceResults.selection.named('entireSurface');

    obj.zetaPlaneResults.label('zetaPlaneResults');
	obj.zetaPlaneResults.selection.named('zetaPlane');
    
    obj.bulkBoundaryResults.label('bulkBoundaryResults');
	obj.bulkBoundaryResults.selection.named('bulkBoundary');
    
    obj.weResults.label('weResults');
	obj.weResults.selection.named('workingElectrode');

    obj.ceResults.label('ceResults');
    obj.ceResults.selection.named('counterElectrode');

    obj.centralCrossectionResults.label('centralCrossectionRseults');
    obj.centralCrossectionResults .set('genpoints', {'0' '0'; '0' '1'});
    
    obj.centralDDLCrossectionResults.label('centralDDLCrossectionResults');
    obj.centralDDLCrossectionResults.set('genpoints', {'0' '0'; '0' num2str(obj.epsilon)});
emd