function createDatasets(obj)
    obj.weResults                   = obj.m.result.dataset.create('weResults', 'Edge2D');
    obj.ceResults                   = obj.m.result.dataset.create('ceResults', 'Edge2D');
    obj.bpeSurfaceResults           = obj.m.result.dataset.create('bpeSurfaceResults', 'Edge2D');
    obj.zetaPlaneResults            = obj.m.result.dataset.create('zetaPlaneResults', 'Edge2D');
    obj.bulkBoundaryResults         = obj.m.result.dataset.create('bulkBoundaryResults', 'Edge2D');    
    obj.entireSurfaceResults        = obj.m.result.dataset.create('entireSurfaceResults', 'Edge2D');
	obj.centralCrossectionResults   = obj.m.result.dataset.create('centralCrossectionResults', 'CutLine2D');
	obj.centralDDLCrossectionResults= obj.m.result.dataset.create('centralDDLCrossectionResults', 'CutLine2D');
    
  	obj.leftBpeEdgeCrossectionResults   = obj.m.result.dataset.create('leftBpeEdgeCrossectionResults', 'CutLine2D');
  	obj.rightBpeEdgeCrossectionResults   = obj.m.result.dataset.create('rightBpeEdgeCrossectionResults', 'CutLine2D');
    
    obj.cathodeCrossectionResults   = obj.m.result.dataset.create('cathodeCrossectionResults', 'CutLine2D');
    obj.anodeCrossectionResults   = obj.m.result.dataset.create('anodeCrossectionResults', 'CutLine2D');

    obj.halfDebyeLengthResults   = obj.m.result.dataset.create('halfDebyeLengthResults', 'CutLine2D');

    obj.m.result.dataset.create('multiPurposeCutLine', 'CutLine2D');
end