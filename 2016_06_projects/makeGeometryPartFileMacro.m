import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

% parameter file for geometry to be created
if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end


caseStudyTitle = 'geometryPartFile';
createEmptyProject;

makeParameterFile
makeVariablesFiles

model = m.m;

%% 1d geometry part

geom_id = 'simple1dGeometryPart';

model.geom.create(geom_id,'Part',1);
model.geom(geom_id).label(geom_id);

model.geom(geom_id).inputParam().set('lambdaD','0.1','Debye length');
% model.geom(geom_id).inputParam().set('L','1','width of electrochemical cell');

model.geom(geom_id).repairTol(1.0E-16);
model.geom(geom_id).create('ddl', 'Interval');
model.geom(geom_id).feature('ddl').set('p','0,lambdaD');
model.geom(geom_id).feature('ddl').set('selresult', 'on');
model.geom(geom_id).feature('ddl').set('selresultshow', 'all');

% model.geom(geom_id).create('space', 'Interval');
% model.geom(geom_id).feature('space').label('space');
% model.geom(geom_id).feature('space').set('p','lambdaD,L');
% model.geom(geom_id).feature('space').set('selresult', 'on');
% model.geom(geom_id).feature('space').set('selresultshow', 'all');

model.geom(geom_id).create('surfaceVertex', 'BoxSelection');
model.geom(geom_id).feature('surfaceVertex').set('entitydim', '0');
model.geom(geom_id).feature('surfaceVertex').label('surfaceVertex');
model.geom(geom_id).feature('surfaceVertex').set('condition', 'inside');
model.geom(geom_id).feature('surfaceVertex').set('xmax', 'lambdaD/2');

model.geom(geom_id).create('bulkVertex', 'BoxSelection');
model.geom(geom_id).feature('bulkVertex').set('entitydim', '0');
model.geom(geom_id).feature('bulkVertex').label('bulkVertex');
model.geom(geom_id).feature('bulkVertex').set('condition', 'inside');
model.geom(geom_id).feature('bulkVertex').set('xmin', 'lambdaD/2');

model.geom(geom_id).create('zetaVertex', 'BoxSelection');
model.geom(geom_id).feature('zetaVertex').set('entitydim', '0');
model.geom(geom_id).feature('zetaVertex').label('zetaVertex');
model.geom(geom_id).feature('zetaVertex').set('condition', 'inside');
model.geom(geom_id).feature('zetaVertex').set('xmin', 'lambdaD/2');
% model.geom(geom_id).feature('zetaVertex').set('xmax', 'L-lambdaD/2');

%% simple rectangular bulk geometry part
geom_id = 'simpleBulkGeometryPart';
% m.updateParameters();
%m.m.geom.create('geom',2);

model.geom.create(geom_id,'Part',2);
model.geom(geom_id).label(geom_id);

model.geom(geom_id).inputParam().set('lambdaD','0.1','Debye length');
model.geom(geom_id).inputParam().set('W','1','width of electrochemical cell');
model.geom(geom_id).inputParam().set('H','1','height of electrochemical cell');
model.geom(geom_id).inputParam().set('XleftBoundary','-0.5','x position of cell''s left boundary');

model.geom(geom_id).localParam().set('XrightBoundary','XleftBoundary+W','x position of cell''s right boundary');

model.geom(geom_id).repairTol(1.0E-10);
model.geom(geom_id).create('space', 'Rectangle');
model.geom(geom_id).feature('space').label('space');
model.geom(geom_id).feature('space').set('selresult', 'on');
model.geom(geom_id).feature('space').set('size', {'W' 'H'});
model.geom(geom_id).feature('space').set('selresultshow', 'all');
model.geom(geom_id).feature('space').set('pos', {'XleftBoundary' '0'});

model.geom(geom_id).create('entireSurface', 'BoxSelection');
model.geom(geom_id).feature('entireSurface').set('entitydim', '1');
model.geom(geom_id).feature('entireSurface').label('entireSurface');
model.geom(geom_id).feature('entireSurface').set('condition', 'inside');
model.geom(geom_id).feature('entireSurface').set('ymax', 'lambdaD/2');
model.geom(geom_id).create('upperBoundary', 'BoxSelection');
model.geom(geom_id).feature('upperBoundary').set('entitydim', '1');
model.geom(geom_id).feature('upperBoundary').label('upperBoundarySelection');
model.geom(geom_id).feature('upperBoundary').set('ymin', 'H-lambdaD/2');
model.geom(geom_id).feature('upperBoundary').set('condition', 'inside');
model.geom(geom_id).create('leftBoundarySelection', 'BoxSelection');
model.geom(geom_id).feature('leftBoundarySelection').set('entitydim', '1');
model.geom(geom_id).feature('leftBoundarySelection').label('leftBoundarySelection');
model.geom(geom_id).feature('leftBoundarySelection').set('condition', 'inside');
model.geom(geom_id).feature('leftBoundarySelection').set('xmax', 'XleftBoundary+lambdaD/2');
model.geom(geom_id).create('rightBoundarySelection', 'BoxSelection');
model.geom(geom_id).feature('rightBoundarySelection').set('entitydim', '1');
model.geom(geom_id).feature('rightBoundarySelection').label('rightBoundarySelection');
model.geom(geom_id).feature('rightBoundarySelection').set('condition', 'inside');
model.geom(geom_id).feature('rightBoundarySelection').set('xmin', 'XrightBoundary-lambdaD/1');
model.geom(geom_id).create('lateralBoundary', 'UnionSelection');
model.geom(geom_id).feature('lateralBoundary').set('entitydim', '1');
model.geom(geom_id).feature('lateralBoundary').label('lateralBoundary');
model.geom(geom_id).feature('lateralBoundary').set('input', {'rightBoundarySelection' 'leftBoundarySelection'});
model.geom(geom_id).create('workingElectrode', 'UnionSelection');
model.geom(geom_id).feature('workingElectrode').set('entitydim', '1');
model.geom(geom_id).feature('workingElectrode').label('workingElectrode');
model.geom(geom_id).feature('workingElectrode').set('input', {'leftBoundarySelection'});
model.geom(geom_id).create('counterElectrode', 'UnionSelection');
model.geom(geom_id).feature('counterElectrode').set('entitydim', '1');
model.geom(geom_id).feature('counterElectrode').label('counterElectrode');
model.geom(geom_id).feature('counterElectrode').set('input', {'rightBoundarySelection'});
model.geom(geom_id).create('electrodes', 'UnionSelection');
model.geom(geom_id).feature('electrodes').set('entitydim', '1');
model.geom(geom_id).feature('electrodes').label('electrodes');
model.geom(geom_id).feature('electrodes').set('input', {'lateralBoundary'});
model.geom(geom_id).create('bulkBoundary', 'UnionSelection');
model.geom(geom_id).feature('bulkBoundary').set('entitydim', '1');
model.geom(geom_id).feature('bulkBoundary').label('bulkBoundary');
model.geom(geom_id).feature('bulkBoundary').set('input', {'lateralBoundary' 'upperBoundary'});
model.geom(geom_id).create('bpeSurface', 'UnionSelection');
model.geom(geom_id).feature('bpeSurface').set('entitydim', '1');
model.geom(geom_id).feature('bpeSurface').label('bpeSurface');
model.geom(geom_id).feature('bpeSurface').set('input', {'entireSurface'});

%     model.geom(geom_id).feature('fin').label('Form Assembly');
%     model.geom(geom_id).feature('fin').set('repairtol', '1.0E-10');
%     model.geom(geom_id).feature('fin').set('action', 'assembly');
%     model.geom(geom_id).run('fin');

%% simple rectangular ddl geometry part
geom_id = 'simpleDdlGeometryPart';
model.geom.create(geom_id,'Part',2);
model.geom(geom_id).label(geom_id);


model.geom(geom_id).inputParam().set('lambdaD','0.1','Debye length');
model.geom(geom_id).inputParam().set('W','1','width of electrochemical cell');
% model.geom(geom_id).inputParam().set('H','1','height of electrochemical cell');
model.geom(geom_id).inputParam().set('XleftBoundary','-0.5','x position of cell''s left boundary');
model.geom(geom_id).localParam().set('XrightBoundary','XleftBoundary+W','x position of cell''s right boundary');


model.geom(geom_id).repairTol(1.0E-10);
model.geom(geom_id).create('ddl', 'Rectangle');
model.geom(geom_id).feature('ddl').label('ddl');
model.geom(geom_id).feature('ddl').set('selresult', 'on');
model.geom(geom_id).feature('ddl').set('size', {'W' 'lambdaD'});
model.geom(geom_id).feature('ddl').set('selresultshow', 'all');
model.geom(geom_id).feature('ddl').set('pos', {'XleftBoundary' '-lambdaD'});

model.geom(geom_id).create('entireSurface', 'BoxSelection');
model.geom(geom_id).feature('entireSurface').set('entitydim', '1');
model.geom(geom_id).feature('entireSurface').label('entireSurface');
model.geom(geom_id).feature('entireSurface').set('condition', 'inside');
model.geom(geom_id).feature('entireSurface').set('ymax', '-lambdaD/2');
model.geom(geom_id).create('upperBoundary', 'BoxSelection');
model.geom(geom_id).feature('upperBoundary').set('entitydim', '1');
model.geom(geom_id).feature('upperBoundary').label('upperBoundarySelection');
model.geom(geom_id).feature('upperBoundary').set('ymin', '-lambdaD/2');
model.geom(geom_id).feature('upperBoundary').set('condition', 'inside');
model.geom(geom_id).create('leftBoundarySelection', 'BoxSelection');
model.geom(geom_id).feature('leftBoundarySelection').set('entitydim', '1');
model.geom(geom_id).feature('leftBoundarySelection').label('leftBoundarySelection');
model.geom(geom_id).feature('leftBoundarySelection').set('condition', 'inside');
model.geom(geom_id).feature('leftBoundarySelection').set('xmax', 'XleftBoundary+lambdaD/2');
model.geom(geom_id).create('rightBoundarySelection', 'BoxSelection');
model.geom(geom_id).feature('rightBoundarySelection').set('entitydim', '1');
model.geom(geom_id).feature('rightBoundarySelection').label('rightBoundarySelection');
model.geom(geom_id).feature('rightBoundarySelection').set('condition', 'inside');
model.geom(geom_id).feature('rightBoundarySelection').set('xmin', 'XrightBoundary-lambdaD/1');
model.geom(geom_id).create('lateralBoundary', 'UnionSelection');
model.geom(geom_id).feature('lateralBoundary').set('entitydim', '1');
model.geom(geom_id).feature('lateralBoundary').label('lateralBoundary');
model.geom(geom_id).feature('lateralBoundary').set('input', {'rightBoundarySelection' 'leftBoundarySelection'});
model.geom(geom_id).create('workingElectrode', 'UnionSelection');
model.geom(geom_id).feature('workingElectrode').set('entitydim', '1');
model.geom(geom_id).feature('workingElectrode').label('workingElectrode');
model.geom(geom_id).feature('workingElectrode').set('input', {'leftBoundarySelection'});
model.geom(geom_id).create('counterElectrode', 'UnionSelection');
model.geom(geom_id).feature('counterElectrode').set('entitydim', '1');
model.geom(geom_id).feature('counterElectrode').label('counterElectrode');
model.geom(geom_id).feature('counterElectrode').set('input', {'rightBoundarySelection'});
model.geom(geom_id).create('electrodes', 'UnionSelection');
model.geom(geom_id).feature('electrodes').set('entitydim', '1');
model.geom(geom_id).feature('electrodes').label('electrodes');
model.geom(geom_id).feature('electrodes').set('input', {'lateralBoundary'});
model.geom(geom_id).create('bulkBoundary', 'UnionSelection');
model.geom(geom_id).feature('bulkBoundary').set('entitydim', '1');
model.geom(geom_id).feature('bulkBoundary').label('bulkBoundary');
model.geom(geom_id).feature('bulkBoundary').set('input', {'lateralBoundary' 'upperBoundary'});
model.geom(geom_id).create('bpeSurface', 'UnionSelection');
model.geom(geom_id).feature('bpeSurface').set('entitydim', '1');
model.geom(geom_id).feature('bpeSurface').label('bpeSurface');
model.geom(geom_id).feature('bpeSurface').set('input', {'entireSurface'});

% model.geom(geom_id).feature('fin').label('Form Assembly');
% model.geom(geom_id).feature('fin').set('repairtol', '1.0E-10');
% model.geom(geom_id).feature('fin').set('action', 'assembly');
% %     model.geom(geom_id).run('fin');

%% geometry sequences
% needs parameters
% makeParameterFile
model.param().loadFile(files('parameterFile'));

% 1d geometry
% bulk geometry
model.modelNode.create('simple1dGeometryComp');
model.geom.create('simple1dGeometry',1);
model.geom('simple1dGeometry').label('simple1dGeometry');
model.geom('simple1dGeometry').repairTol(1.0E-16);
model.geom('simple1dGeometry').create('simple1dGeometryPartInstance','PartInstance');
model.geom('simple1dGeometry').feature('simple1dGeometryPartInstance').set('part','simple1dGeometryPart');
model.geom('simple1dGeometry').feature('simple1dGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
% model.geom('simple1dGeometry').feature('simple1dGeometryPartInstance').setEntry('inputexpr', 'L', '(1+extendedDdlFactor)*lambdaD');
model.geom('simple1dGeometry').feature('fin').label('Form Union');
model.geom('simple1dGeometry').feature('fin').set('repairtol', '1.0E-16');
model.geom('simple1dGeometry').feature('fin').set('action', 'union');
model.geom('simple1dGeometry').run('fin');

% bulk geometry
model.modelNode.create('simpleBulkGeometryComp');
model.geom.create('simpleBulkGeometry',2);
model.geom('simpleBulkGeometry').label('simpleBulkGeometry');
model.geom('simpleBulkGeometry').repairTol(1.0E-10);
model.geom('simpleBulkGeometry').create('simpleBulkGeometryPartInstance','PartInstance');
model.geom('simpleBulkGeometry').feature('simpleBulkGeometryPartInstance').set('part','simpleBulkGeometryPart');
model.geom('simpleBulkGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
model.geom('simpleBulkGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'W', 'W');
model.geom('simpleBulkGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'H', 'H');
model.geom('simpleBulkGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'XleftBoundary', 'XleftBoundary');
model.geom('simpleBulkGeometry').feature('fin').label('Form Assembly');
model.geom('simpleBulkGeometry').feature('fin').set('repairtol', '1.0E-10');
model.geom('simpleBulkGeometry').feature('fin').set('action', 'assembly');
model.geom('simpleBulkGeometry').run('fin');

% ddl geometry
model.modelNode.create('simpleDdlGeometryComp');
model.geom.create('simpleDdlGeometry',2);
model.geom('simpleDdlGeometry').label('simpleDdlGeometry');
model.geom('simpleDdlGeometry').repairTol(1.0E-10);
model.geom('simpleDdlGeometry').create('simpleDdlGeometryPartInstance','PartInstance');
model.geom('simpleDdlGeometry').feature('simpleDdlGeometryPartInstance').set('part','simpleDdlGeometryPart');
model.geom('simpleDdlGeometry').feature('simpleDdlGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
model.geom('simpleDdlGeometry').feature('simpleDdlGeometryPartInstance').setEntry('inputexpr', 'W', 'W');
model.geom('simpleDdlGeometry').feature('simpleDdlGeometryPartInstance').setEntry('inputexpr', 'XleftBoundary', 'XleftBoundary');
model.geom('simpleDdlGeometry').feature('fin').label('Form Assembly');
model.geom('simpleDdlGeometry').feature('fin').set('repairtol', '1.0E-10');
model.geom('simpleDdlGeometry').feature('fin').set('action', 'assembly');
model.geom('simpleDdlGeometry').run('fin');

% simple assembled geometry
model.modelNode.create('simpleAssembledGeometryComp');
model.geom.create('simpleAssembledGeometry',2);
model.geom('simpleAssembledGeometry').label('simpleAssembledGeometry');
model.geom('simpleAssembledGeometry').repairTol(1.0E-10);
model.geom('simpleAssembledGeometry').create('simpleBulkGeometryPartInstance','PartInstance');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').label('simpleBulkGeometryPartInstance');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').set('part','simpleBulkGeometryPart');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'W', 'W');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'H', 'H');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'XleftBoundary', 'XleftBoundary');
model.geom('simpleAssembledGeometry').create('simpleDdlGeometryPartInstance','PartInstance');
model.geom('simpleAssembledGeometry').feature('simpleDdlGeometryPartInstance').label('simpleDdlGeometryPartInstance');
model.geom('simpleAssembledGeometry').feature('simpleDdlGeometryPartInstance').set('part','simpleDdlGeometryPart');
model.geom('simpleAssembledGeometry').feature('simpleDdlGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
model.geom('simpleAssembledGeometry').feature('simpleDdlGeometryPartInstance').setEntry('inputexpr', 'W', 'W');
model.geom('simpleAssembledGeometry').feature('simpleDdlGeometryPartInstance').setEntry('inputexpr', 'XleftBoundary', 'XleftBoundary');
model.geom('simpleAssembledGeometry').create('entireSurface', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('entireSurface').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('entireSurface').label('entireSurface');
model.geom('simpleAssembledGeometry').feature('entireSurface').set('input', {'simpleDdlGeometryPartInstance_entireSurface'});
model.geom('simpleAssembledGeometry').create('upperBoundary', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('upperBoundary').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('upperBoundary').label('upperBoundarySelection');
model.geom('simpleAssembledGeometry').feature('upperBoundary').set('input', {'simpleBulkGeometryPartInstance_upperBoundary'});
model.geom('simpleAssembledGeometry').create('leftBoundarySelection', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('leftBoundarySelection').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('leftBoundarySelection').label('leftBoundarySelection');
model.geom('simpleAssembledGeometry').feature('leftBoundarySelection').set('input', {'simpleBulkGeometryPartInstance_leftBoundarySelection', 'simpleDdlGeometryPartInstance_leftBoundarySelection'});
model.geom('simpleAssembledGeometry').create('rightBoundarySelection', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('rightBoundarySelection').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('rightBoundarySelection').label('rightBoundarySelection');
model.geom('simpleAssembledGeometry').feature('rightBoundarySelection').set('input', {'simpleBulkGeometryPartInstance_rightBoundarySelection', 'simpleDdlGeometryPartInstance_rightBoundarySelection'});
model.geom('simpleAssembledGeometry').create('lateralBoundary', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('lateralBoundary').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('lateralBoundary').label('lateralBoundary');
model.geom('simpleAssembledGeometry').feature('lateralBoundary').set('input', {'rightBoundarySelection' 'leftBoundarySelection'});
model.geom('simpleAssembledGeometry').create('workingElectrode', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('workingElectrode').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('workingElectrode').label('workingElectrode');
model.geom('simpleAssembledGeometry').feature('workingElectrode').set('input', {'leftBoundarySelection'});
model.geom('simpleAssembledGeometry').create('counterElectrode', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('counterElectrode').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('counterElectrode').label('counterElectrode');
model.geom('simpleAssembledGeometry').feature('counterElectrode').set('input', {'rightBoundarySelection'});
model.geom('simpleAssembledGeometry').create('electrodes', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('electrodes').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('electrodes').label('electrodes');
model.geom('simpleAssembledGeometry').feature('electrodes').set('input', {'lateralBoundary'});
model.geom('simpleAssembledGeometry').create('bulkBoundary', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('bulkBoundary').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('bulkBoundary').label('bulkBoundary');
model.geom('simpleAssembledGeometry').feature('bulkBoundary').set('input', {'lateralBoundary' 'upperBoundary'});
model.geom('simpleAssembledGeometry').create('bpeSurface', 'UnionSelection');
model.geom('simpleAssembledGeometry').feature('bpeSurface').set('entitydim', '1');
model.geom('simpleAssembledGeometry').feature('bpeSurface').label('bpeSurface');
model.geom('simpleAssembledGeometry').feature('bpeSurface').set('input', {'entireSurface'});

model.geom('simpleAssembledGeometry').feature('fin').label('Form Assembly');
model.geom('simpleAssembledGeometry').feature('fin').set('repairtol', '1.0E-10');
model.geom('simpleAssembledGeometry').feature('fin').set('action', 'assembly');
model.geom('simpleAssembledGeometry').run('fin');

%% geometries with mesh refinement

% bulk geometry
model.modelNode.create('simpleBulkGeometryWithMeshingEntitiesComp');
model.geom.create('simpleBulkGeometryWithMeshingEntities',2);
model.geom('simpleBulkGeometryWithMeshingEntities').label('simpleBulkGeometryWithMeshingEntities');
model.geom('simpleBulkGeometryWithMeshingEntities').repairTol(1.0E-10);
model.geom('simpleBulkGeometryWithMeshingEntities').create('simpleBulkGeometryPartInstance','PartInstance');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').set('part','simpleBulkGeometryPart');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'W', 'W');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'H', 'H');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'XleftBoundary', 'XleftBoundary');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('fin').label('Form Assembly');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('fin').set('repairtol', '1.0E-10');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('fin').set('action', 'assembly');
% model.geom('simpleBulkGeometryWithMeshingEntities').run('fin');

% model.geom('simpleBulkGeometryWithMeshingEntities').selection.create('csel1', 'CumulativeSelection');
% model.geom('simpleBulkGeometryWithMeshingEntities').selection('csel1').label('leftBpeBoundaryEdge');
% model.geom('simpleBulkGeometryWithMeshingEntities').selection.create('csel2', 'CumulativeSelection');
% model.geom('simpleBulkGeometryWithMeshingEntities').selection('csel2').label('rightBpeBoundaryEdge');
% model.geom('simpleBulkGeometryWithMeshingEntities').create('simpleBulkGeometryPartInstance', 'PartInstance');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').set('inputexpr', {'lambdaD' 'W' 'H' 'XleftBoundary'});
model.geom('simpleBulkGeometryWithMeshingEntities').feature('simpleBulkGeometryPartInstance').set('selkeepnoncontr', true);
model.geom('simpleBulkGeometryWithMeshingEntities').create('leftBpeBoundaryEdge', 'Polygon');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('leftBpeBoundaryEdge').set('selresult', 'on');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('leftBpeBoundaryEdge').set('contributeto', 'csel1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('leftBpeBoundaryEdge').set('x', '-Wbpe/2,-Wbpe/2');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('leftBpeBoundaryEdge').set('y', '0,H');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('leftBpeBoundaryEdge').set('type', 'open');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('leftBpeBoundaryEdge').set('selresultshow', 'bnd');
model.geom('simpleBulkGeometryWithMeshingEntities').create('rightBpeBoundaryEdge', 'Polygon');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('rightBpeBoundaryEdge').set('selresult', 'on');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('rightBpeBoundaryEdge').set('contributeto', 'csel2');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('rightBpeBoundaryEdge').set('x', 'Wbpe/2,Wbpe/2');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('rightBpeBoundaryEdge').set('y', '0,H');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('rightBpeBoundaryEdge').set('type', 'open');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('rightBpeBoundaryEdge').set('selresultshow', 'bnd');
model.geom('simpleBulkGeometryWithMeshingEntities').create('pt1', 'Point');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt1').setIndex('p', '-Wbpe/2-L', 0, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt1').setIndex('p', '0', 1, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt1').set('selresult', 'on');
model.geom('simpleBulkGeometryWithMeshingEntities').create('pt2', 'Point');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt2').setIndex('p', 'Wbpe/2-L', 0, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt2').setIndex('p', '0', 1, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt2').set('selresult', 'on');
model.geom('simpleBulkGeometryWithMeshingEntities').create('pt3', 'Point');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt3').setIndex('p', '-Wbpe/2+L', 0, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt3').setIndex('p', '0', 1, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt3').set('selresult', 'on');
model.geom('simpleBulkGeometryWithMeshingEntities').create('pt4', 'Point');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt4').setIndex('p', '+Wbpe/2+L', 0, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt4').setIndex('p', '0', 1, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt4').set('selresult', 'on');
model.geom('simpleBulkGeometryWithMeshingEntities').create('pt5', 'Point');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt5').setIndex('p', '-Wbpe/2', 0, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt5').setIndex('p', 'L', 1, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt5').set('selresult', 'on');
model.geom('simpleBulkGeometryWithMeshingEntities').create('pt6', 'Point');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt6').setIndex('p', 'Wbpe/2', 0, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt6').setIndex('p', 'L', 1, 0);
model.geom('simpleBulkGeometryWithMeshingEntities').feature('pt6').set('selresult', 'on');
% 
% model.geom('simpleBulkGeometryWithMeshingEntities').create('unisel1', 'UnionSelection');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('unisel1').set('entitydim', '0');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('unisel1').label('surfacePartitioner');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('unisel1').set('input', {'pt4' 'pt3' 'pt2' 'pt1'});

model.geom('simpleBulkGeometryWithMeshingEntities').create('meshControlEdges', 'UnionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('meshControlEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('meshControlEdges').label('meshControlEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('meshControlEdges').set('input', {'leftBpeBoundaryEdge' 'rightBpeBoundaryEdge'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('westwardMeshControlVertices', 'UnionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('westwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('westwardMeshControlVertices').label('westwardMeshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('westwardMeshControlVertices').set('input', {'pt2' 'pt1'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('eastwardMeshControlVertices', 'UnionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('eastwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('eastwardMeshControlVertices').label('eastwardMeshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('eastwardMeshControlVertices').set('input', {'pt4' 'pt3'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('meshControlVertices', 'UnionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('meshControlVertices').set('entitydim', '0');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('meshControlVertices').label('meshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('meshControlVertices').set('input', {'pt6' 'pt5' 'pt4' 'pt3' 'pt2' 'pt1'});

model.geom('simpleBulkGeometryWithMeshingEntities').feature('fin').set('repairtol', '1.0E-10');
model.geom('simpleBulkGeometryWithMeshingEntities').run('fin');

model.geom('simpleBulkGeometryWithMeshingEntities').create('boundariesAdjacentToMeshControlVertices', 'AdjacentSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToMeshControlVertices').set('entitydim', '0');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToMeshControlVertices').label('boundariesAdjacentToMeshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToMeshControlVertices').set('input', {'meshControlVertices'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundariesAdjacentToMeshControlEdges', 'AdjacentSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToMeshControlEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToMeshControlEdges').label('boundariesAdjacentToMeshControlEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToMeshControlEdges').set('input', {'meshControlEdges'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundariesAdjacentToWestwardMeshControlVertices', 'AdjacentSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToWestwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToWestwardMeshControlVertices').label('boundariesAdjacentToWestwardMeshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToWestwardMeshControlVertices').set('input', {'westwardMeshControlVertices'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundariesAdjacentToEastwardMeshControlVertices', 'AdjacentSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToEastwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToEastwardMeshControlVertices').label('boundariesAdjacentToEastwardMeshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToEastwardMeshControlVertices').set('input', {'eastwardMeshControlVertices'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundaryAtBpeRegion', 'IntersectionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryAtBpeRegion').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryAtBpeRegion').label('boundaryAtBpeRegion');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryAtBpeRegion').set('input', {'boundariesAdjacentToEastwardMeshControlVertices' 'boundariesAdjacentToWestwardMeshControlVertices'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('westwardThinningEdges', 'IntersectionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('westwardThinningEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('westwardThinningEdges').label('westwardThinningEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('westwardThinningEdges').set('input', {'boundariesAdjacentToEastwardMeshControlVertices' 'boundariesAdjacentToMeshControlEdges'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('eastwardThinningEdges', 'IntersectionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('eastwardThinningEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('eastwardThinningEdges').label('eastwardThinningEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('eastwardThinningEdges').set('input', {'boundariesAdjacentToMeshControlEdges' 'boundariesAdjacentToWestwardMeshControlVertices'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('horizontalThinningEdges', 'UnionSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('horizontalThinningEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('horizontalThinningEdges').label('horizontalThinningEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('horizontalThinningEdges').set('input', {'eastwardThinningEdges' 'westwardThinningEdges'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundariesAdjacentToHorizontalMeshControlVertices', 'AdjacentSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToHorizontalMeshControlVertices').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToHorizontalMeshControlVertices').label('boundariesAdjacentToHorizontalMeshControlVertices');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToHorizontalMeshControlVertices').set('input', {'horizontalThinningEdges'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('lateralThinningEdges', 'DifferenceSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('lateralThinningEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('lateralThinningEdges').label('lateralThinningEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('lateralThinningEdges').set('add', {'boundariesAdjacentToHorizontalMeshControlVertices'});
model.geom('simpleBulkGeometryWithMeshingEntities').feature('lateralThinningEdges').set('subtract', {'simpleBulkGeometryPartInstance_entireSurface'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundariesAdjacentToLateralThinningEdges', 'AdjacentSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToLateralThinningEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToLateralThinningEdges').label('boundariesAdjacentToLateralThinningEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundariesAdjacentToLateralThinningEdges').set('input', {'lateralThinningEdges'});
model.geom('simpleBulkGeometryWithMeshingEntities').create('boundaryOfBpePartitioningEdges', 'DifferenceSelection');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryOfBpePartitioningEdges').set('entitydim', '1');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryOfBpePartitioningEdges').label('boundaryOfBpePartitioningEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryOfBpePartitioningEdges').set('add', {'boundariesAdjacentToLateralThinningEdges'});
model.geom('simpleBulkGeometryWithMeshingEntities').feature('boundaryOfBpePartitioningEdges').set('subtract', {'horizontalThinningEdges'});
% model.geom('simpleBulkGeometryWithMeshingEntities').create('mcv1', 'MeshControlVertices');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('mcv1').active(false);
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('mcv1').selection('input').named('meshControlVertices');
% model.geom('simpleBulkGeometryWithMeshingEntities').create('mce1', 'MeshControlEdges');
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('mce1').active(false);
% model.geom('simpleBulkGeometryWithMeshingEntities').feature('mce1').selection('input').named('meshControlEdges');
model.geom('simpleBulkGeometryWithMeshingEntities').run;

% ddl geometry
model.modelNode.create('simpleDdlGeometryWithMeshingEntitiesComp');
model.geom.create('simpleDdlGeometryWithMeshingEntities',2);
model.geom('simpleDdlGeometryWithMeshingEntities').label('simpleDdlGeometryWithMeshingEntities');
model.geom('simpleDdlGeometryWithMeshingEntities').repairTol(1.0E-10);
model.geom('simpleDdlGeometryWithMeshingEntities').create('simpleDdlGeometryPartInstance', 'PartInstance');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('simpleDdlGeometryPartInstance').label('simpleDdlGeometryPartInstance');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('simpleDdlGeometryPartInstance').set('part', 'simpleDdlGeometryPart');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('simpleDdlGeometryPartInstance').set('inputexpr', {'lambdaD' 'W' 'XleftBoundary'});
model.geom('simpleDdlGeometryWithMeshingEntities').feature('simpleDdlGeometryPartInstance').set('selkeepnoncontr', true);
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt1', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt1').setIndex('p', '-Wbpe/2-L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt1').setIndex('p', '0', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt1').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt2', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt2').setIndex('p', '-Wbpe/2-L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt2').setIndex('p', '-lambdaD', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt2').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt3', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt3').setIndex('p', 'Wbpe/2-L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt3').setIndex('p', '0', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt3').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt4', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt4').setIndex('p', 'Wbpe/2-L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt4').setIndex('p', '-lambdaD', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt4').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt5', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt5').setIndex('p', '-Wbpe/2+L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt5').setIndex('p', '0', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt5').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt6', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt6').setIndex('p', '-Wbpe/2+L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt6').setIndex('p', '-lambdaD', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt6').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt7', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt7').setIndex('p', '+Wbpe/2+L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt7').setIndex('p', '0', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt7').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt8', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt8').setIndex('p', '+Wbpe/2+L', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt8').setIndex('p', '-lambdaD', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt8').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt9', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt9').setIndex('p', '-Wbpe/2', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt9').setIndex('p', '0', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt9').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt10', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt10').setIndex('p', '-Wbpe/2', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt10').setIndex('p', '-lambdaD', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt10').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt11', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt11').setIndex('p', 'Wbpe/2', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt11').setIndex('p', '0', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt11').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pt12', 'Point');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt12').setIndex('p', 'Wbpe/2', 0, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt12').setIndex('p', '-lambdaD', 1, 0);
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pt12').set('selresult', 'on');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pol1', 'Polygon');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol1').set('x', '-Wbpe/2-L,-Wbpe/2-L');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol1').set('y', '-lambdaD,0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol1').set('type', 'open');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pol2', 'Polygon');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol2').set('x', 'Wbpe/2-L,Wbpe/2-L');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol2').set('y', '-lambdaD,0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol2').set('type', 'open');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pol3', 'Polygon');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol3').set('x', '-Wbpe/2+L,-Wbpe/2+L');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol3').set('y', '-lambdaD,0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol3').set('type', 'open');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pol4', 'Polygon');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol4').set('x', 'Wbpe/2+L,Wbpe/2+L');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol4').set('y', '-lambdaD,0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol4').set('type', 'open');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pol5', 'Polygon');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol5').set('x', '-Wbpe/2,-Wbpe/2');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol5').set('y', '-lambdaD,0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol5').set('type', 'open');
model.geom('simpleDdlGeometryWithMeshingEntities').create('pol6', 'Polygon');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol6').set('x', 'Wbpe/2,Wbpe/2');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol6').set('y', '-lambdaD,0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('pol6').set('type', 'open');
model.geom('simpleDdlGeometryWithMeshingEntities').create('westwardMeshControlVertices', 'UnionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('westwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('westwardMeshControlVertices').label('westwardMeshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('westwardMeshControlVertices').set('input', {'pt4' 'pt3' 'pt2' 'pt1'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('eastwardMeshControlVertices', 'UnionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('eastwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('eastwardMeshControlVertices').label('eastwardMeshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('eastwardMeshControlVertices').set('input', {'pt8' 'pt7' 'pt6' 'pt5'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('centralMeshControlVertices', 'UnionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('centralMeshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('centralMeshControlVertices').label('centralMeshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('centralMeshControlVertices').set('input', {'pt12' 'pt11' 'pt10' 'pt9'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('meshControlVertices', 'UnionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('meshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('meshControlVertices').label('meshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('meshControlVertices').set('input', {'centralMeshControlVertices' 'eastwardMeshControlVertices' 'westwardMeshControlVertices'});

model.geom('simpleDdlGeometryWithMeshingEntities').feature('fin').set('repairtol', '1.0E-10');
model.geom('simpleDdlGeometryWithMeshingEntities').run('fin');

model.geom('simpleDdlGeometryWithMeshingEntities').create('boundariesAdjacentToCentralMeshControlVertices', 'AdjacentSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToCentralMeshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToCentralMeshControlVertices').label('boundariesAdjacentToCentralMeshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToCentralMeshControlVertices').set('input', {'centralMeshControlVertices'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('boundariesAdjacentToWestwardMeshControlVertices', 'AdjacentSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToWestwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToWestwardMeshControlVertices').label('boundariesAdjacentToWestwardMeshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToWestwardMeshControlVertices').set('input', {'westwardMeshControlVertices'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('boundariesAdjacentToEastwardMeshControlVertices', 'AdjacentSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToEastwardMeshControlVertices').set('entitydim', '0');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToEastwardMeshControlVertices').label('boundariesAdjacentToEastwardMeshControlVertices');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('boundariesAdjacentToEastwardMeshControlVertices').set('input', {'eastwardMeshControlVertices'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('eastwardThinningEdges', 'IntersectionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('eastwardThinningEdges').set('entitydim', '1');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('eastwardThinningEdges').label('eastwardThinningEdges');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('eastwardThinningEdges').set('input', {'boundariesAdjacentToWestwardMeshControlVertices' 'boundariesAdjacentToCentralMeshControlVertices'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('westwardThinningEdges', 'IntersectionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('westwardThinningEdges').set('entitydim', '1');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('westwardThinningEdges').label('westwardThinningEdges');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('westwardThinningEdges').set('input', {'boundariesAdjacentToCentralMeshControlVertices' 'boundariesAdjacentToEastwardMeshControlVertices'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('horizontalThinningEdges', 'UnionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('horizontalThinningEdges').set('entitydim', '1');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('horizontalThinningEdges').label('horizontalThinningEdges');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('horizontalThinningEdges').set('input', {'westwardThinningEdges' 'eastwardThinningEdges'});
model.geom('simpleDdlGeometryWithMeshingEntities').create('lateralThinningEdges', 'BoxSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('lateralThinningEdges').set('entitydim', '1');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('lateralThinningEdges').label('lateralThinningEdges');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('lateralThinningEdges').set('ymin', '-lambdaD/2');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('lateralThinningEdges').set('ymax', '-lambdaD/2');
model.geom('simpleDdlGeometryWithMeshingEntities').create('centralBpeMeshControlEdges', 'IntersectionSelection');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('centralBpeMeshControlEdges').set('entitydim', '1');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('centralBpeMeshControlEdges').label('centralBpeMeshControlEdges');
model.geom('simpleDdlGeometryWithMeshingEntities').feature('centralBpeMeshControlEdges').set('input', {'boundariesAdjacentToEastwardMeshControlVertices' 'boundariesAdjacentToWestwardMeshControlVertices'});

%% meshing
if exist('hMaxFactor','var')
    m.hMaxFactor = hMaxFactor;
else
    m.hMaxFactor = 0.1;
end

m.mesh1D();

a0 = m.intFirstDebyeLength(1);
r = m.intFirstDebyeLength(2)/m.intFirstDebyeLength(1);
N = ceil(log(1-m.epsilon/a0*(1-r))/log(r));
R = max(r^N, m.intFirstDebyeLength(end)/m.intFirstDebyeLength(1));

a_lat = m.intFirstDebyeLength(end);
% N_lat = ceil(m.w_mesh / a_lat);

% N_lat_end = ceil((obj.w_mesh-obj.epsilon) / a_lat);

% jumped 2 debye lengths
% hmax = m.intRemaining(end);
% 
% r_zeta = m.intRemaining(1)/m.intFirstDebyeLength(end);
% a_zeta = m.intRemaining(1);
% 
% r_bulk = m.intRemaining(1)/m.intFirstDebyeLength(end);
% a0_bulk = m.intRemaining(1);
% aN_bulk = a0_bulk*R^N;

hmax = m.intRemaining(end);

r_zeta = m.intFirstDebyeLength(end)/m.intFirstDebyeLength(end-1);
a_zeta = m.intExtendedDdl(1);

r_bulk = m.intRemaining(2)/m.intRemaining(1);
a0_bulk = m.intRemaining(1);
aN_bulk = a0_bulk*r_bulk^N;

% 1d geometry 
model.mesh.create('simple1dGeometryRefinedMesh', 'simple1dGeometry');
model.mesh('simple1dGeometryRefinedMesh').create('edg1', 'Edge'); % ddl
model.mesh('simple1dGeometryRefinedMesh').feature('edg1').create('dis1', 'Distribution');
model.mesh('simple1dGeometryRefinedMesh').feature('edg1').create('size1', 'Size');

model.mesh('simple1dGeometryRefinedMesh').feature('edg1').selection.geom(1);
model.mesh('simple1dGeometryRefinedMesh').feature('edg1').selection.all;

model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('dis1').selection.named('simple1dGeometry_simple1dGeometryPartInstance_ddl_dom');
model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('dis1').set('type','explicit');
model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('dis1').set('explicit', m.distributionFirstDebyeLengthStr);

% model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('size1').selection.named('simple1dGeometry_simple1dGeometryPartInstance_zetaVertex');
% model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('size1').set('custom', 'on');
% model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('size1').set('hmaxactive', true);
% model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('size1').set('hgrad', r_zeta);
% model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('size1').set('hmax', sprintf('%e*L',a_zeta));
% model.mesh('simple1dGeometryRefinedMesh').feature('edg1').feature('size1').set('hgradactive', true);

model.mesh('simple1dGeometryRefinedMesh').feature('size').set('custom', 'on');
model.mesh('simple1dGeometryRefinedMesh').feature('size').set('hmin', sprintf('%e*L',a0));
model.mesh('simple1dGeometryRefinedMesh').feature('size').set('hgrad', r_zeta);
model.mesh('simple1dGeometryRefinedMesh').feature('size').set('hmax', sprintf('%e*L',a0_bulk));

model.mesh('simple1dGeometryRefinedMesh').run;
model.mesh('simple1dGeometryRefinedMesh').export.set('type', 'nativeascii');
simple1dGeometryRefinedMeshFile = [pwd,'\',m.projectPath,'\simple1dGeometryRefinedMesh.mphtxt'];
files('simple1dGeometryRefinedMeshFile') = simple1dGeometryRefinedMeshFile;
model.mesh('simple1dGeometryRefinedMesh').export(simple1dGeometryRefinedMeshFile);

% bulk geometry 
model.mesh.create('simpleBulkGeometryRefinedMesh', 'simpleBulkGeometryWithMeshingEntities');
model.mesh('simpleBulkGeometryRefinedMesh').create('edg1', 'Edge'); % horizontal thinning edges
model.mesh('simpleBulkGeometryRefinedMesh').create('edg2', 'Edge'); % lateral
model.mesh('simpleBulkGeometryRefinedMesh').create('ftri1', 'FreeTri');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').create('dis1', 'Distribution');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').create('dis2', 'Distribution');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg2').create('dis1', 'Distribution');
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').create('size1', 'Size'); % bpe surface refinement

model.mesh('simpleBulkGeometryRefinedMesh').label('bulkGeometryPartRefinedMesh');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').selection.named('simpleBulkGeometryWithMeshingEntities_horizontalThinningEdges');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis1').selection.named('simpleBulkGeometryWithMeshingEntities_eastwardThinningEdges');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis2').selection.named('simpleBulkGeometryWithMeshingEntities_westwardThinningEdges');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg2').selection.named('simpleBulkGeometryWithMeshingEntities_lateralThinningEdges');
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').selection.named('simpleBulkGeometryWithMeshingEntities_boundaryAtBpeRegion'); %
% model.mesh('simpleBulkGeometryRefinedMesh').feature('size').set('hauto', 3);
model.mesh('simpleBulkGeometryRefinedMesh').feature('size').set('custom', 'on');
% model.mesh('simpleBulkGeometryRefinedMesh').feature('size').set('hgrad', '1.1');
model.mesh('simpleBulkGeometryRefinedMesh').feature('size').set('hgrad', r_bulk);
model.mesh('simpleBulkGeometryRefinedMesh').feature('size').set('hmin',  sprintf('%e*L',a0_bulk));
model.mesh('simpleBulkGeometryRefinedMesh').feature('size').set('hmax', 'H/4');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis1').set('elemcount', N);
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis1').set('method', 'geometric');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis1').set('type', 'predefined');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis1').set('elemratio', R);
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis2').set('elemcount', N);
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis2').set('method', 'geometric');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis2').set('reverse', true);
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis2').set('type', 'predefined');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg1').feature('dis2').set('elemratio', R);
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg2').feature('dis1').set('elemcount', N);
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg2').feature('dis1').set('type', 'predefined');
model.mesh('simpleBulkGeometryRefinedMesh').feature('edg2').feature('dis1').set('elemratio', R);
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('custom', 'on');
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('hmaxactive', true);
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('hgrad', r_bulk);
% model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('hmin',  sprintf('%e*L',a_zeta)); % minimum element width
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('hmax', 'H/4'); % maximum element width at bpe surface
model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('hgradactive', true);
% model.mesh('simpleBulkGeometryRefinedMesh').feature('ftri1').feature('size1').set('hminactive', false);
model.mesh('simpleBulkGeometryRefinedMesh').run;

model.mesh('simpleBulkGeometryRefinedMesh').export.set('type', 'nativeascii');
simpleBulkGeometryRefinedMeshFile = [pwd,'\',m.projectPath,'\simpleBulkGeometryRefinedMesh.mphtxt'];
% model.mesh('simpleBulkGeometryRefinedMesh').export.set('filename', 'D:\windows\Users\jotelha\Google Docs\johnny\matlab\jlhbpe\2016_06_projects\dat\2016_06_13_18_52_56_geometryPartFile\simpeDdlGeometryRefinedMesh.mphtxt');
files('simpleBulkGeometryRefinedMeshFile') = simpleBulkGeometryRefinedMeshFile;
model.mesh('simpleBulkGeometryRefinedMesh').export(simpleBulkGeometryRefinedMeshFile);

% ddl geometry

model.mesh.create('simpleDdlGeometryRefinedMesh', 'simpleDdlGeometryWithMeshingEntities');
model.mesh('simpleDdlGeometryRefinedMesh').create('edg1', 'Edge'); % horizontal thinning edges
model.mesh('simpleDdlGeometryRefinedMesh').create('edg2', 'Edge'); % lateral thinning edges
model.mesh('simpleDdlGeometryRefinedMesh').create('edg3', 'Edge'); % centralBpeMeshControl, seems to be empty, dont know why
model.mesh('simpleDdlGeometryRefinedMesh').create('edg4', 'Edge'); % remaining
model.mesh('simpleDdlGeometryRefinedMesh').create('map1', 'Map');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').create('dis1', 'Distribution');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').create('dis2', 'Distribution');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').create('dis1', 'Distribution');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').create('size2', 'Size');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').create('size1', 'Size');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').create('size1', 'Size');

model.mesh('simpleDdlGeometryRefinedMesh').label('simpleDdlGeometryRefinedMesh');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').selection.named('simpleDdlGeometryWithMeshingEntities_horizontalThinningEdges');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis1').selection.named('simpleDdlGeometryWithMeshingEntities_eastwardThinningEdges');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis2').selection.named('simpleDdlGeometryWithMeshingEntities_westwardThinningEdges');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').selection.named('simpleDdlGeometryWithMeshingEntities_lateralThinningEdges');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').feature('dis1').selection.named('simpleDdlGeometryWithMeshingEntities_lateralThinningEdges');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').selection.named('simpleDdlGeometryWithMeshingEntities_centralBpeMeshControlEdges');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size1').selection.named('simpleDdlGeometryWithMeshingEntities_centralBpeMeshControlEdges');
% model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size2').selection.geom('simpleDdlGeometryWithMeshingEntities', 1);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size2').selection.named('simpleDdlGeometryWithMeshingEntities_meshControlVertices');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').selection.remaining;
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').feature('size1').selection.named('simpleDdlGeometryWithMeshingEntities_meshControlVertices');
model.mesh('simpleDdlGeometryRefinedMesh').feature('size').set('custom', 'on');
model.mesh('simpleDdlGeometryRefinedMesh').feature('size').set('hmin', sprintf('%e*L',a0));
model.mesh('simpleDdlGeometryRefinedMesh').feature('size').set('hgrad', r_zeta);
model.mesh('simpleDdlGeometryRefinedMesh').feature('size').set('hmax', 'H/4');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis1').set('elemcount', N);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis1').set('method', 'geometric');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis1').set('type', 'predefined');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis1').set('elemratio', R);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis2').set('elemcount', N);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis2').set('method', 'geometric');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis2').set('reverse', true);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis2').set('type', 'predefined');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg1').feature('dis2').set('elemratio', R);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').feature('dis1').set('elemcount', N);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').feature('dis1').set('method', 'geometric');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').feature('dis1').set('type', 'predefined');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg2').feature('dis1').set('elemratio', R);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size2').set('custom', 'on');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size2').set('hmaxactive', true);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size2').set('hmax', 'L');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size1').set('custom', 'on');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size1').set('hmaxactive', true);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size1').set('hgrad', r_zeta);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size1').set('hmax', 'H/4');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg3').feature('size1').set('hgradactive', true);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').feature('size1').set('custom', 'on');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').feature('size1').set('hmaxactive', true);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').feature('size1').set('hgrad', r_zeta);
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').feature('size1').set('hmax', 'L');
model.mesh('simpleDdlGeometryRefinedMesh').feature('edg4').feature('size1').set('hgradactive', true);
model.mesh('simpleDdlGeometryRefinedMesh').run;

model.mesh('simpleDdlGeometryRefinedMesh').export.set('type', 'nativeascii');
simpleDdlGeometryRefinedMeshFile = [pwd,'\',m.projectPath,'\simpleDdlGeometryRefinedMesh.mphtxt'];
files('simpleDdlGeometryRefinedMeshFile') = simpleDdlGeometryRefinedMeshFile;
model.mesh('simpleDdlGeometryRefinedMesh').export(simpleDdlGeometryRefinedMeshFile);

m.saveState;

%% save information on created files
geometryPartsProjectName = m.projectName;
geometryPartsMphFile = char(m.m.getFilePath);
files('geometryPartsMphFile') = geometryPartsMphFile;

jlh.hf.saveMapAsTxt(files,'globalFiles.txt');
% toSave = {'geometryPartsProjectName','geometryPartsMphFile','simpleBulkGeometryRefinedMeshFile','simpleDdlGeometryRefinedMeshFile'};
% N = numel(toSave);
% dat = cell(N,1);
% for i=1:N
%     dat{i} = eval(toSave{i});
% end
% T = cell2table(dat,'RowNames',toSave);
% writetable(T,'geometryPartsFile.txt','WriteRowNames',true,'WriteVariableNames',false,'Delimiter',' ');
% save('geometryPartsFile.mat',toSave{:},'-mat');

% simpleDdlGeometryProjectName = m.projectName;
% simpleDdlGeometryMphFile = m.m.getFilePath;

% ModelUtil.remove(m.model_tag);