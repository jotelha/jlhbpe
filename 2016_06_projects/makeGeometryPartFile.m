import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

% parameter file for geometry to be created
caseStudyParameterFile = 'parameters_duval2001bipolar.m';


caseStudyTitle = 'geometryPartFile';
createEmptyProject;

% simple rectangular bulk geometry part
geom_id = 'simpleBulkGeometryPart';
% m.updateParameters();
%m.m.geom.create('geom',2);

model = m.m;
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
m.updateParameters;

model.modelNode.create('comp');

% bulk geometry
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
model.geom.create('simpleAssembledGeometry',2);
model.geom('simpleAssembledGeometry').label('simpleAssembledGeometry');
model.geom('simpleAssembledGeometry').repairTol(1.0E-10);
model.geom('simpleAssembledGeometry').create('simpleBulkGeometryPartInstance','PartInstance');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').set('part','simpleBulkGeometryPart');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'lambdaD', 'lambdaD');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'W', 'W');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'H', 'H');
model.geom('simpleAssembledGeometry').feature('simpleBulkGeometryPartInstance').setEntry('inputexpr', 'XleftBoundary', 'XleftBoundary');
model.geom('simpleAssembledGeometry').create('simpleDdlGeometryPartInstance','PartInstance');
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
model.geom('simpleAssembledGeometry').feature('entireSurface').set('input', {'simpleBulkGeometryPartInstance_entireSurface'});
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

m.saveState;

%% save information on created files
geometryPartsProjectName = m.projectName;
geometryPartsMphFile = char(m.m.getFilePath);

toSave = {'geometryPartsProjectName','geometryPartsMphFile'};
N = numel(toSave);
dat = cell(N,1);
for i=1:N
    dat{i} = eval(toSave{i});
end
T = cell2table(dat,'RowNames',toSave,'VariableNames',{''});
writetable(T,'geometryPartsFile.txt','WriteRowNames',true,'WriteVariableNames',false);
save('geometryPartsFile.mat',toSave{:},'-mat');

% simpleDdlGeometryProjectName = m.projectName;
% simpleDdlGeometryMphFile = m.m.getFilePath;

% ModelUtil.remove(m.model_tag);
