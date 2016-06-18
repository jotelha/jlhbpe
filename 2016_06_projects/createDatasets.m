% model.result.dataset.create('weResults', 'Edge2D');
% model.result.dataset.create('ceResults', 'Edge2D');
model.result.dataset.create('bpeSurfaceResults', 'Edge2D');
model.result.dataset.create('zetaPlaneResults', 'Edge2D');
model.result.dataset.create('bulkBoundaryResults', 'Edge2D');    
model.result.dataset.create('entireSurfaceResults', 'Edge2D');
model.result.dataset.create('centralCrossectionResults', 'CutLine2D');
model.result.dataset.create('centralDDLCrossectionResults', 'CutLine2D');
    
model.result.dataset.create('leftBpeEdgeCrossectionResults', 'CutLine2D');
model.result.dataset.create('rightBpeEdgeCrossectionResults', 'CutLine2D');
    
model.result.dataset.create('cathodeCrossectionResults', 'CutLine2D');
model.result.dataset.create('anodeCrossectionResults', 'CutLine2D');

model.result.dataset.create('halfDebyeLengthResults', 'CutLine2D');

model.result.dataset.create('multiPurposeCutLine', 'CutLine2D');