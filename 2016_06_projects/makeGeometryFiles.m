import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

% parameter file for geometry to be created
caseStudyParameterFile = 'parameters_duval2001bipolar.m';


%% simple bulk geometry
caseStudyTitle = 'simpleBulkGeometry';
createEmptyProject;

% m.m.modelNode.create('comp'); 
m.updateParameters();
%m.m.geom.create('geom',2);
m.m.geom.create('geom','Part',2);

    model = m.m;
    model.geom(geom_id).repairTol(1.0E-10);
%     model.geom(geom_id).create('ddl', 'Rectangle');
%     model.geom(geom_id).feature('ddl').active(false);
%     model.geom(geom_id).feature('ddl').label('ddl');
%     model.geom(geom_id).feature('ddl').set('selresult', 'on');
%     model.geom(geom_id).feature('ddl').set('size', {'w' 'epsilon'});
%     model.geom(geom_id).feature('ddl').set('selresultshow', 'all');
%     model.geom(geom_id).feature('ddl').set('pos', {'-w_bpe/2-w_bulkLeft' '0'});
    model.geom(geom_id).create('space', 'Rectangle');
    model.geom(geom_id).feature('space').label('space');
    model.geom(geom_id).feature('space').set('selresult', 'on');
    model.geom(geom_id).feature('space').set('size', {'W' 'L'});
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
    model.geom(geom_id).create('workingElectrode', 'UnionSelection');
    model.geom(geom_id).feature('workingElectrode').set('entitydim', '1');
    model.geom(geom_id).feature('workingElectrode').label('workingElectrode');
    model.geom(geom_id).create('counterElectrode', 'UnionSelection');
    model.geom(geom_id).feature('counterElectrode').set('entitydim', '1');
    model.geom(geom_id).feature('counterElectrode').label('counterElectrode');
    model.geom(geom_id).create('electrodes', 'UnionSelection');
    model.geom(geom_id).feature('electrodes').set('entitydim', '1');
    model.geom(geom_id).feature('electrodes').label('electrodes');
    model.geom(geom_id).create('bulkBoundary', 'UnionSelection');
    model.geom(geom_id).feature('bulkBoundary').set('entitydim', '1');
    model.geom(geom_id).feature('bulkBoundary').label('bulkBoundary');
    model.geom(geom_id).create('bpeSurface', 'UnionSelection');
    model.geom(geom_id).feature('bpeSurface').set('entitydim', '1');
    model.geom(geom_id).feature('bpeSurface').label('bpeSurface');
    model.geom(geom_id).feature('bpeSurface').set('input', {'entireSurface'});
    
    model.geom(geom_id).feature('fin').label('Form Assembly');
    model.geom(geom_id).feature('fin').set('repairtol', '1.0E-10');
    model.geom(geom_id).feature('fin').set('action', 'assembly');
%     model.geom(geom_id).run('fin');

m.saveState;

simpleBulkGeometryProjectName = m.projectName;
simpleBulkGeometryMphFile = m.m.getFilePath;

ModelUtil.remove(m.model_tag);


%% simple ddl geometry
caseStudyTitle = 'simpleDdlGeometry';
createEmptyProject;
% m.m.modelNode.create('comp'); 
m.updateParameters();
% m.m.geom.create('geom',2);
m.m.geom.create('geom','Part',2);
% m.makeSimpleRectangularDdlGeometry('geom');

    model = m.m;
    model.geom(geom_id).repairTol(1.0E-10);

    model.geom(geom_id).create('ddl', 'Rectangle');
    model.geom(geom_id).feature('ddl').label('ddl');
    model.geom(geom_id).feature('ddl').set('selresult', 'on');
    model.geom(geom_id).feature('ddl').set('size', {'W' 'lambdaD'});
    model.geom(geom_id).feature('ddl').set('selresultshow', 'all');
    model.geom(geom_id).feature('ddl').set('pos', {'XleftBoundary' '-lambdaD'});

%     model.geom(geom_id).create('space', 'Rectangle');
%     model.geom(geom_id).feature('space').label('space');
%     model.geom(geom_id).feature('space').set('selresult', 'on');
%     model.geom(geom_id).feature('space').set('size', {'W' 'L'});
%     model.geom(geom_id).feature('space').set('selresultshow', 'all');
%     model.geom(geom_id).feature('space').set('pos', {'XleftBoundary' '0'});
   
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
    model.geom(geom_id).create('workingElectrode', 'UnionSelection');
    model.geom(geom_id).feature('workingElectrode').set('entitydim', '1');
    model.geom(geom_id).feature('workingElectrode').label('workingElectrode');
    model.geom(geom_id).create('counterElectrode', 'UnionSelection');
    model.geom(geom_id).feature('counterElectrode').set('entitydim', '1');
    model.geom(geom_id).feature('counterElectrode').label('counterElectrode');
    model.geom(geom_id).create('electrodes', 'UnionSelection');
    model.geom(geom_id).feature('electrodes').set('entitydim', '1');
    model.geom(geom_id).feature('electrodes').label('electrodes');
    model.geom(geom_id).create('bulkBoundary', 'UnionSelection');
    model.geom(geom_id).feature('bulkBoundary').set('entitydim', '1');
    model.geom(geom_id).feature('bulkBoundary').label('bulkBoundary');
    model.geom(geom_id).create('bpeSurface', 'UnionSelection');
    model.geom(geom_id).feature('bpeSurface').set('entitydim', '1');
    model.geom(geom_id).feature('bpeSurface').label('bpeSurface');
    model.geom(geom_id).feature('bpeSurface').set('input', {'entireSurface'});
    
    model.geom(geom_id).feature('fin').label('Form Assembly');
    model.geom(geom_id).feature('fin').set('repairtol', '1.0E-10');
    model.geom(geom_id).feature('fin').set('action', 'assembly');
%     model.geom(geom_id).run('fin');

m.saveState;

simpleDdlGeometryProjectName = m.projectName;
simpleDdlGeometryMphFile = m.m.getFilePath;

ModelUtil.remove(m.model_tag);

%% simple assembled geometry
caseStudyTitle = 'simpleAssembledGeometry';
createEmptyProject;
m.m.modelNode.create('comp'); 
m.updateParameters();
m.m.geom.create('geom',2);

m.m.geom('geom').repairTol(1.0E-10);
m.m.geom('geom').insertFile(simpleBulkGeometryMphFile,'geom');
m.m.geom('geom').insertFile(simpleDdlGeometryMphFile,'geom');

m.m.geom('geom').feature('fin').label('Form Assembly');
m.m.geom('geom').feature('fin').set('repairtol', '1.0E-10');
m.m.geom('geom').feature('fin').set('action', 'assembly');

m.saveState;

simpleAssembledGeometryProjectName = m.projectName;
simpleAssembledGeometryMphFile = m.m.getFilePath;

ModelUtil.remove(m.model_tag);