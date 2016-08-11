import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end

caseStudyTitle = 'dses2dDirichlet';
createEmptyProject;

% or, to load an existing project from server or from file
% tag = 'xyz'
% loadExistingProject

% prepare text files
% makeParameterFile
% makeVariablesFiles

if ~exist('BatchDir','var')
%     BatchDir = [ pwd, '\', m.projectPath, '\batch'];
    BatchDir = 'batch';
end
if ~exist('ParametricDir','var')
%     BatchDir = [ pwd, '\', m.projectPath, '\batch'];
    ParametricDir = 'parametric';
end
if ~exist('OutputDir','var')
%     BatchDir = [ pwd, '\', m.projectPath, '\batch'];
    OutputDir = 'output';
end
 
BatchDirAbsolute = [ pwd(), '\', m.projectPath, '\', BatchDir];
ParametricDirAbsolute = [pwd(), '\', m.projectPath, '\', ParametricDir];
OutputDirAbsolute = [pwd(),'\',m.projectPath, '\', OutputDir];

if ~exist(BatchDirAbsolute,'dir')
    mkdir(BatchDirAbsolute);
end
if ~exist(ParametricDirAbsolute,'dir')
    mkdir(ParametricDirAbsolute);
end
if ~exist(OutputDirAbsolute,'dir')
    mkdir(OutputDirAbsolute);
end

GlobalEvaluationsTableFileName = [OutputDir,'\globalEvaluations.txt'];
exportTertiaryCurrentDistributionDataFile = [OutputDir,'\exportTertiaryCurrentDistributionDataFile.txt'];
exportTertiaryCurrentDistribution2dDataFile = [OutputDir,'\exportTertiaryCurrentDistribution2dDataFile.txt'];

%% load parameters, geometry and mesh parts
% T = readtable('geometryPartsFile.txt','ReadRowNames',true,'ReadVariableNames',false,'Delimiter',' ');
T = readtable('globalFiles.txt','ReadRowNames',false,'ReadVariableNames',false,'Delimiter',' ');
% writetable(T,'geometryPartsFile.txt','WriteRowNames',true,'WriteVariableNames',false,'Delimiter',' ');
C = table2cell(T);
files = containers.Map(C(:,1),C(:,2));
% geometryPartsMphFile = T.(1)('geometryPartsMphFile');
% simpleBulkGeometryRefinedMeshFile = T.(1)('simpleBulkGeometryRefinedMeshFile');
% simpleDdlGeometryRefinedMeshFile = T.(1)('simpleDdlGeometryRefinedMeshFile');

model = m.m;

model.param().loadFile(files('parameterFile'));

%% mesh components
model.modelNode.create('mcomp1', 'MeshComponent');
model.geom.create('mgeom1', 2); % mesh geometry
model.mesh.create('simpleBulkGeometryRefinedMeshPart', 'mgeom1');
model.mesh('simpleBulkGeometryRefinedMeshPart').create('imp1', 'Import');
model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('filename', files('simpleBulkGeometryRefinedMeshFile'));
model.mesh('simpleBulkGeometryRefinedMeshPart').run;
% model.mesh('mpart1').feature('imp1').set('facepartition', 'auto');

model.modelNode.create('mcomp2', 'MeshComponent');
model.geom.create('mgeom2', 2); % mesh geometry
model.mesh.create('simpleDdlGeometryRefinedMeshPart', 'mgeom2');
model.mesh('simpleDdlGeometryRefinedMeshPart').create('imp1', 'Import');
model.mesh('simpleDdlGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
model.mesh('simpleDdlGeometryRefinedMeshPart').feature('imp1').set('filename', files('simpleDdlGeometryRefinedMeshFile'));
model.mesh('simpleDdlGeometryRefinedMeshPart').run;

%% create
% component for rough tertiary current approximation
% m.comp_id = 'tertiaryCurrentDistributionComponent';
model.modelNode.create('dilutedSpeciesAndElectrostatics2dComponent'); 
model.geom.create('dilutedSpeciesAndElectrostatics2dGeometry',2);
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').insertFile(files('geometryPartsMphFile'), 'simpleAssembledGeometry');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('fin').set('repairtol', '1e-10');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('fin').set('action', 'assembly');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('fin').set('createpairs', 'off');

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').create('workingElectrode', 'UnionSelection');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').create('counterElectrode', 'UnionSelection');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').create('electrodes', 'UnionSelection');

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('workingElectrode').set('entitydim', '1');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('workingElectrode').label('workingElectrodeUnion');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('workingElectrode').set('input', {'simpleDdlGeometryPartInstance1_workingElectrode', 'simpleBulkGeometryPartInstance1_workingElectrode'});

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('counterElectrode').set('entitydim', '1');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('counterElectrode').label('counterElectrodeUnion');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('counterElectrode').set('input', {'simpleDdlGeometryPartInstance1_counterElectrode', 'simpleBulkGeometryPartInstance1_counterElectrode'});

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('electrodes').set('entitydim', '1');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('electrodes').label('electrodesUnion');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('electrodes').set('input', {'workingElectrode', 'counterElectrode'});

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').create('ddlOuterBoundary', 'UnionSelection');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').create('bulkOuterBoundary', 'UnionSelection');

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('ddlOuterBoundary').set('entitydim', '1');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('ddlOuterBoundary').label('ddlOuterBoundary');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('ddlOuterBoundary').set('input', {'simpleDdlGeometryPartInstance1_electrodes', 'simpleDdlGeometryPartInstance1_entireSurface'});

model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('bulkOuterBoundary').set('entitydim', '1');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('bulkOuterBoundary').label('bulkOuterBoundary');
model.geom('dilutedSpeciesAndElectrostatics2dGeometry').feature('bulkOuterBoundary').set('input', {'simpleBulkGeometryPartInstance1_electrodes', 'simpleBulkGeometryPartInstance1_upperBoundary'});

% pair
model.pair.create('p1', 'Identity', 'dilutedSpeciesAndElectrostatics2dGeometry');

% mesh
model.mesh.create('dilutedSpeciesAndElectrostatics2dMesh', 'dilutedSpeciesAndElectrostatics2dGeometry');
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').create('copy1', 'Copy');
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').create('copy2', 'Copy');

% functions
model.func.create('smoothenBpeBC', 'Rectangle');
model.func.create('interpolateStoredValues1d','Interpolation');
model.func.create('interpolateStoredValues2d','Interpolation');


% analytical pb
% model.func.create('phi_pb', 'Analytic');
% model.func.create('phi_pbx', 'Analytic');
% for i = 1:m.numberOfSpecies   
%     model.func.create(m.c_pb_id{i}, 'Analytic');
% end
% model.func.create('pbDecayFunction', 'Analytic');

%operators
model.cpl.create('extrudeZetaPlaneToSurface', 'LinearExtrusion', 'dilutedSpeciesAndElectrostatics2dGeometry');

% model.cpl.create('onSurface', 'LinearExtrusion', 'dilutedSpeciesAndElectrostatics2dGeometry');
% model.cpl.create('intSurface', 'Integration', 'dilutedSpeciesAndElectrostatics2dGeometry');
% model.cpl.create('intBulk', 'Integration', 'dilutedSpeciesAndElectrostatics2dGeometry');

% variables
model.variable.create('domainVariables');
model.variable.create('surfaceVariables');
model.variable.create('weVariables');
model.variable.create('ceVariables');
% model.variable.create('bulkVariables');

% physics
model.physics.create('DilutedSpecies', 'DilutedSpecies', 'dilutedSpeciesAndElectrostatics2dGeometry', m.c_id);
model.physics.create('Electrostatics', 'Electrostatics', 'dilutedSpeciesAndElectrostatics2dGeometry');

model.physics('DilutedSpecies').feature.create('init2', 'init', 2);
model.physics('Electrostatics').feature.create('init2', 'init', 2);

model.physics('DilutedSpecies').feature.create('BulkConcentration', 'Concentration', 1);
model.physics('DilutedSpecies').feature.create('DdlConcentration', 'Concentration', 1);
% model.physics('TertiaryCurrentDistribution').feature.create('BpeSurface', 'ExternalElectrodeSurface', 1);
% model.physics('DilutedSpecies').feature.create('SurfaceFlux', 'Fluxes', 1);
% model.physics('DilutedSpecies').feature.create('BulkFlux', 'Fluxes', 0);

model.physics('Electrostatics').feature.create('SpaceChargeDensity', 'SpaceChargeDensity', 2);
% model.physics('Electrostatics').feature.create('SurfaceChargeDensity', 'SurfaceChargeDensity', 1);
% model.physics('Electrostatics').feature.create('BulkPotential', 'ElectricPotential', 1);
model.physics('Electrostatics').feature.create('BulkPotential', 'ElectricPotential', 1);
model.physics('Electrostatics').feature.create('DdlPotential', 'ElectricPotential', 1);

model.physics('DilutedSpecies').create('DilutedSpeciesContinuity', 'Continuity', 1);
model.physics('Electrostatics').create('ElectrostaticsContinuity', 'Continuity', 1);

% for i=1:m.nReactions
% %     model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature.create(m.reactionNames{i}, 'ElectrodeReaction', 1);
% end

%% update

% dimensional space charge density term
chargeDensitySummand = prepTerm('z_id*c_id','z_id','c_id',m.z_id,m.c_id);
chargeDensityTerm = strcat('F_const*(',... % 
        strjoin(chargeDensitySummand','+' ), ')');
    
pbDecayFunction = 'exp(-(x/lambdaD+delta))';

pdePSourceTerm = prepTerm('chargeDensityRampFactor*chargeDensityTerm','chargeDensityTerm',chargeDensityTerm); % allow for ramping

% pair
model.pair('p1').source.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_upperBoundary');
model.pair('p1').destination.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleBulkGeometryPartInstance1_entireSurface');

% mesh
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy1').selection('source').geom(2);
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy1').selection('destination').geom(2);
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy1').selection('source').all;
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy1').selection('destination').named('dilutedSpeciesAndElectrostatics2dGeometry_simpleBulkGeometryPartInstance1_space_dom');
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').run('copy1');

model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy2').set('mesh', 'simpleDdlGeometryRefinedMeshPart');
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy2').selection('source').geom(2);
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy2').selection('destination').geom(2);
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy2').selection('source').all;
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').feature('copy2').selection('destination').named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_ddl_dom');
model.mesh('dilutedSpeciesAndElectrostatics2dMesh').run('copy2');

% functions 
model.func('smoothenBpeBC').set('upper', 'w_bpe/2');
model.func('smoothenBpeBC').set('smooth', 'epsilon*smootheningFactor');
model.func('smoothenBpeBC').set('funcname', 'smoothenBpeBC');
model.func('smoothenBpeBC').set('lower', '-w_bpe/2');

% interpolate previous ddl results
model.func('interpolateStoredValues1d').set('source', 'file');
model.func('interpolateStoredValues1d').set('interp', 'linear');

model.func('interpolateStoredValues1d').setIndex('funcs', 'phi_interp_ddl', 0, 0);
model.func('interpolateStoredValues1d').setIndex('funcs', num2str(1), 0, 1);

for i=1:m.numberOfSpecies
    model.func('interpolateStoredValues1d').setIndex('funcs', sprintf('%s_interp_ddl',m.c_id{i}), i, 0);
    model.func('interpolateStoredValues1d').setIndex('funcs', num2str(i+1), i, 1);
end
% model.func('interpolateStoredValues1d').set('filename', files('exportDilutedSpeciesAndElectrostatics1dDataFile'));
% model.func('interpolateStoredValues1d').set('nargs', '2');
% model.func('interpolateStoredValues1d').importData;

% interpolate previous bulk results
model.func('interpolateStoredValues2d').set('source', 'file');
model.func('interpolateStoredValues2d').set('interp', 'linear');

model.func('interpolateStoredValues2d').setIndex('funcs', 'phi_interp_bulk', 0, 0);
model.func('interpolateStoredValues2d').setIndex('funcs', num2str(1), 0, 1);

for i=1:m.numberOfSpecies
    model.func('interpolateStoredValues2d').setIndex('funcs', sprintf('%s_interp_bulk',m.c_id{i}), i, 0);
    model.func('interpolateStoredValues2d').setIndex('funcs', num2str(i+1), i, 1);
end
% model.func('interpolateStoredValues2d').set('filename', files('exportTertiaryCurrentDistribution2dDataFile'));
% model.func('interpolateStoredValues2d').set('nargs', '2');
% model.func('interpolateStoredValues2d').importData;


% model.func('pbDecayFunction').set('funcname', 'pbDecayFunction');
% model.func('pbDecayFunction').set('args', {'x'});
% model.func('pbDecayFunction').set('expr', pbDecayFunction);
% model.func('pbDecayFunction').set('fununit', '1');
% model.func('pbDecayFunction').set('argunit', 'm');
% model.func('pbDecayFunction').set('plotargs', {'x' '0' 'L'});

% operators
model.cpl('extrudeZetaPlaneToSurface').set('opname','extrudeZetaPlaneToSurface');
% model.cpl('onSurface').selection('srcvertex2').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
% model.cpl('onSurface').selection('srcvertex1').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
% model.cpl('onSurface').selection('dstvertex2').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
% model.cpl('onSurface').selection('dstvertex1').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
model.cpl('extrudeZetaPlaneToSurface').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 1);
model.cpl('extrudeZetaPlaneToSurface').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_upperBoundary');
model.cpl('extrudeZetaPlaneToSurface').selection('srcvertex1').set(2); % order of vertex creation important
model.cpl('extrudeZetaPlaneToSurface').selection('srcvertex2').set(4);
model.cpl('extrudeZetaPlaneToSurface').selection('dstvertex1').set(1);
model.cpl('extrudeZetaPlaneToSurface').selection('dstvertex2').set(3);
% model.cpl('onSurface').set('opname','onSurface');
% model.cpl('onSurface').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 0);
% model.cpl('onSurface').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
% model.cpl('onSurface').selection('srcvertex1').set(1); % order of vertex creation important
% model.cpl('onSurface').selection('dstvertex1').set(3);

% model.cpl('intSurface').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 0);
% model.cpl('intSurface').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
% model.cpl('intSurface').set('opname', 'intSurface');

% model.cpl('intBulk').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 0);
% model.cpl('intBulk').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simple1dGeometryPartInstance1_bulkVertex');
% model.cpl('intBulk').set('opname', 'intBulk');

% variables
model.variable('domainVariables').model('dilutedSpeciesAndElectrostatics2dComponent');
model.variable('domainVariables').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 2);
model.variable('domainVariables').selection.all;
model.variable('domainVariables').loadFile(files('domainVariablesDilutedSpeciesAndElectrostatics2d'));

model.variable('surfaceVariables').model('dilutedSpeciesAndElectrostatics2dComponent');
model.variable('surfaceVariables').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 1);
model.variable('surfaceVariables').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_bpeSurface');
model.variable('surfaceVariables').loadFile(files('surfaceVariablesDilutedSpeciesAndElectrostatics2d'));

model.variable('weVariables').model('dilutedSpeciesAndElectrostatics2dComponent');
model.variable('weVariables').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 1);
model.variable('weVariables').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_workingElectrode');
model.variable('weVariables').loadFile(files('weVariablesDilutedSpeciesAndElectrostatics2d'));

model.variable('ceVariables').model('dilutedSpeciesAndElectrostatics2dComponent');
model.variable('ceVariables').selection.geom('dilutedSpeciesAndElectrostatics2dGeometry', 1);
model.variable('ceVariables').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_counterElectrode');
model.variable('ceVariables').loadFile(files('ceVariablesDilutedSpeciesAndElectrostatics2d'));

% physics

model.physics('DilutedSpecies').prop('TransportMechanism').set('Convection', false);
model.physics('DilutedSpecies').prop('TransportMechanism').set('Migration', true);

model.physics('Electrostatics').field('electricpotential').field('phi');
model.physics('DilutedSpecies').field('concentration').component(m.c_id);

model.physics('Electrostatics').prop('ShapeProperty').set('order_electricpotential', '5');
model.physics('Electrostatics').prop('ShapeProperty').set('valueType', 'real');

model.physics('DilutedSpecies').prop('ShapeProperty').set('order_concentration', '5');
model.physics('DilutedSpecies').prop('ShapeProperty').set('valueType', 'real');

% standard settings
model.physics('DilutedSpecies').feature('cdm1').set('minput_temperature', 'T');
model.physics('DilutedSpecies').feature('cdm1').set('V', 'phi');
model.physics('Electrostatics').feature('ccn1').set('epsilonr_mat', 'userdef');
model.physics('Electrostatics').feature('ccn1').set('epsilonr', {'epsilon_r' '0' '0' '0' 'epsilon_r' '0' '0' '0' 'epsilon_r'});

model.physics('DilutedSpecies').feature('init2').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_ddl_dom');
model.physics('Electrostatics').feature('init2').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_ddl_dom');

model.physics('Electrostatics').feature('SpaceChargeDensity').selection.all;
model.physics('Electrostatics').feature('SpaceChargeDensity').set('rhoq', pdePSourceTerm); % space charge density

% potential initial values
model.physics('Electrostatics').feature('init1').set('phi', 'phi_interp_bulk(x,y)');
model.physics('Electrostatics').feature('init2').set('phi', 'phi_interp_ddl(lambdaD+y,x)');

% bc selection
% model.physics('DilutedSpecies').feature('SurfaceFlux').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_bpeSurface');
% model.physics('DilutedSpecies').feature('BulkFlux').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simple1dGeometryPartInstance1_bulkVertex');
% model.physics('Electrostatics').feature('SurfaceChargeDensity').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_bpeSurface');

model.physics('DilutedSpecies').feature('BulkConcentration').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_bulkOuterBoundary');
model.physics('DilutedSpecies').feature('DdlConcentration').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_ddlOuterBoundary');
model.physics('Electrostatics').feature('BulkPotential').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_bulkOuterBoundary');
model.physics('Electrostatics').feature('DdlPotential').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_ddlOuterBoundary');
% model.physics('Electrostatics').feature('BulkPotential').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simple1dGeometryPartInstance1_bulkVertex');
% model.physics('Electrostatics').feature('ElectrodePotential').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_electrodes');

% model.physics('Electrostatics').feature('SurfaceChargeDensity').set('rhoqs', 'smoothenBpeBC(X/L)*C_Stern_dimensional*(phi_s-phi)');
model.physics('Electrostatics').feature('BulkPotential').set('V0', 'phi_interp_bulk(x,y)');
model.physics('Electrostatics').feature('DdlPotential').set('V0', 'phi_interp_ddl(lambdaD+y,x)'); % feeder electrodes

for i = 1:m.numberOfSpecies
    D_c_id = sprintf('D_%s',m.c_id{i});
    % isotropic diffusivity
    model.physics('DilutedSpecies').feature('cdm1').set(D_c_id, {m.D_id{i} '0' '0' '0' m.D_id{i} '0' '0' '0' m.D_id{i}});

    model.physics('DilutedSpecies').feature('cdm1').setIndex('z', m.z_id{i}, i-1);

    % bulk and initial concentrations
    
    model.physics('DilutedSpecies').feature('BulkConcentration').setIndex('species', true, i-1);
    model.physics('DilutedSpecies').feature('BulkConcentration').setIndex('c0', sprintf('%s_interp_bulk(x,y)',m.c_id{i}), i-1);
    model.physics('DilutedSpecies').feature('DdlConcentration').setIndex('species', true, i-1);
    model.physics('DilutedSpecies').feature('DdlConcentration').setIndex('c0', sprintf('%s_interp_ddl(lambdaD+y,x)',m.c_id{i}), i-1);
    model.physics('DilutedSpecies').feature('init1').setIndex('initc', sprintf('%s_interp_bulk(x,y)',m.c_id{i}), i-1);
    model.physics('DilutedSpecies').feature('init2').setIndex('initc', sprintf('%s_interp_ddl(lambdaD+y,x)',m.c_id{i}), i-1);

    
%     model.physics('DilutedSpecies').feature('SurfaceFlux').setIndex('species', true, i-1);
%     model.physics('DilutedSpecies').feature('SurfaceFlux').setIndex('N0', sprintf('smoothenBpeBC(X/L)*%s', m.N_id{i}), i-1);
%  
%     model.physics('DilutedSpecies').feature('BulkFlux').setIndex('species', true, i-1);
%     model.physics('DilutedSpecies').feature('BulkFlux').setIndex('N0', sprintf('-smoothenBpeBC(X/L)*onSurface(%s)', m.N_id{i}), i-1);
%  
end

model.physics('DilutedSpecies').feature('DilutedSpeciesContinuity').setIndex('pairs', 'p1', 0);
model.physics('Electrostatics').feature('ElectrostaticsContinuity').setIndex('pairs', 'p1', 0);

%% make batch study
model.result.table.create('tbl1', 'Table'); % needs to be created ahead

if ~exist('pname','var') || ~exist('plistarr','var')
    pname = { 'sweep' 'DeltaPHI' };
    plistarr = { '1,2,3,4' '1.5,2.5,3.5,4.5' };
end

% ParametricSolutionFileName = [BatchDir, '\', m.model_tag, '_parametric.mph']; % absolute path!
ParametricSolutionFileName = 'parametric.mph'; % absolute path!
ParametricSolutionFileNameAbsolute = [ParametricDirAbsolute , '\', ParametricSolutionFileName]; % absolute path!

% BatchFileName = [m.model_tag, '_batch.mph']; % only file name, no path
BatchFileName = 'batch.mph'; % only file name, no path
BatchSolutionFileNameAbsolute = [BatchDirAbsolute , '\', BatchFileName]; % absolute path!

% pname = {'phi_bpe_init','phi_bpe','surfaceFluxRampFactor', 'deltaPhi'};
% plistarr = cellstr( cellfun(@(c) num2str(c), sweepParameterValues,'UniformOutput',false));
if ~exist('pname','var')
    pname = {'DeltaPHI','surfaceFluxRampFactor'};
% pval = cellstr( cellfun(@(c) num2str(c), sweepPhiInitValues,'UniformOutput',false));
    plistarr = {'1'; '1'};
end

punit = '';
useparam = 'on';
pcontinuationmode = 'no';
          
studyNum = model.study.size() + 1;
study_id = sprintf('study%u',studyNum);

currentStationaryStudy = model.study.create(study_id);

solNum = model.sol.size() + 1;
sol_id = sprintf('sol%u',solNum);

model.sol.create(sol_id);
model.sol(sol_id).study(study_id); % ?
model.sol(sol_id).attach(study_id); % ?

% for i=1:nStudySteps
i = 1;
fprintf('\t Creating step %d and setting standard parameters...\n',i);
studyStep_id = sprintf('stationaryStudyStep%d',i);
compile_id = sprintf('st%d',i);
variables_id = sprintf('v%d',i);
solver_id = sprintf('s%d',i);
store_id = sprintf('su%d',i);
batchStep_id = sprintf('batch%d',i);
paramStep_id = sprintf('param%d',i);
batch_id = sprintf('b%d',i);
param_id = sprintf('p%d',i);

model.study(study_id).create(studyStep_id, 'Stationary');

model.sol(sol_id).create(compile_id, 'StudyStep');
model.sol(sol_id).create(variables_id, 'Variables');
model.sol(sol_id).create(solver_id, 'Stationary');
model.sol(sol_id).create(store_id, 'StoreSolution');

% parametric 
model.study(study_id).create(paramStep_id,'Parametric');
model.batch.create(param_id,'Parametric');

model.batch(param_id).create('cl1','Class'); % exterior java class to call for adaption before solving
model.batch(param_id).create('cl2','Class'); % exterior java class to call for adaption before solving
model.batch(param_id).create('so1', 'Solutionseq');
model.batch(param_id).create('saDef', 'Save');
% model.batch(param_id).create('en1', 'Evalnumericalseq');
model.batch(param_id).create('nu1', 'Numericalseq');
model.batch(param_id).create('ex1', 'Exportseq');
model.batch(param_id).study(study_id);


model.study(study_id).feature(paramStep_id).set('pname', pname);
model.study(study_id).feature(paramStep_id).set('paramselect', 'off');
model.study(study_id).feature(paramStep_id).set('plistarr', plistarr);
model.study(study_id).feature(paramStep_id).set('punit', {''});
model.study(study_id).feature(paramStep_id).set('save', 'on'); % important for saving solutions in seperate files
model.study(study_id).feature(paramStep_id).set('keepsol', 'last'); % important for saving solutions in seperate files
model.study(study_id).feature(paramStep_id).set('filename', [ParametricDir, '\', ParametricSolutionFileName]);

model.batch(param_id).set('punit', {''});
model.batch(param_id).set('err', true);
model.batch(param_id).set('keeplog', true);
model.batch(param_id).set('plistarr', plistarr);
model.batch(param_id).set('pname', pname);
model.batch(param_id).set('control', paramStep_id);
model.batch(param_id).feature('cl1').set('filename', 'UpdateInitialValues.class');
model.batch(param_id).feature('cl1').set('input', {'sweep' 'input\inputFiles1d.txt' 'input' 'interpolateStoredValues1d'});
model.batch(param_id).feature('cl2').set('filename', 'UpdateInitialValues.class');
model.batch(param_id).feature('cl2').set('input', {'sweep' 'input\inputFiles2d.txt' 'input' 'interpolateStoredValues2d'});
model.batch(param_id).feature('so1').set('seq', sol_id);
model.batch(param_id).feature('nu1').set('table', 'tbl1');

model.batch(param_id).attach(study_id);

% batch
model.study(study_id).create(batchStep_id, 'Batch');
model.batch.create(batch_id,'Batch');

model.batch(batch_id).create('jo1', 'Jobseq');
% model.batch(batch_id).feature('daDef').create('pr1', 'Process');
model.batch(batch_id).study(study_id);

model.study(study_id).feature(batchStep_id).set('batchdir', BatchDir);
model.study(study_id).feature(batchStep_id).set('batchfile', BatchFileName);
model.study(study_id).feature(batchStep_id).set('specserverdir', 'off');


model.batch(batch_id).set('control', batchStep_id);
model.batch(batch_id).set('batchfile', BatchFileName);
model.batch(batch_id).set('batchdir', BatchDir);
model.batch(batch_id).set('speccomsoldir', 'off');

% model.batch('b1').feature('daDef').feature('pr1').set('clientfilename', 'batch\batch_test_study_batch.mph');
% model.batch('b1').feature('daDef').feature('pr1').set('cmd', {'comsolbatch' '-study' 'std1' '-alivetime' '15' '-inputfile' 'batch_test.mph' '-outputfile' '".\batch\batch_test_study_batch.mph"' '-batchlog'  ...
% '".\batch\batch_test_study_batch.mph.log"'});
% model.batch('b1').feature('daDef').feature('pr1').set('filename', 'batch\batch_test_study_batch.mph');

model.batch(batch_id).attach(study_id);

% set scaling for variables
model.sol(sol_id).feature(variables_id).set('scalemethod','init');

model.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
% model.study(study_id).feature(studyStep_id).set('showdistribute', true);

% model.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);
model.sol(sol_id).feature(solver_id).set('probesel', 'none');
model.sol(sol_id).feature(solver_id).set('nonlin', 'on');
model.sol(sol_id).feature(solver_id).set('stol', '1e-3');

% model.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'const');
model.sol(sol_id).feature(solver_id).feature('dDef').set('ooc', 'on');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
% model.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
% model.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermconst', 'itertol');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '15');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
%% numerical evaluations
% model.result.table('tbl1').comments('Global Evaluation 1 (intWE(tcdee.Ilx))');
model.result.export.create('tbl1', 'Table');
model.result.export('tbl1').set('filename', GlobalEvaluationsTableFileName);

titles = {'DeltaPHI','PHI_bpe', 'Ix_WE', 'Ix_CE', 'Iy_BPE', ' I_total', 'I_anodic', 'I_cathodic', 'I_faradaic', 'I_ohmic'};
expressions = { 'DeltaPHI', ...
                'intWE(Ix)', ...
                'intCE(Ix)', ...
                'intBPE(Iy)', ...
                'intBPE(i_total)',...
                'intBPE(i_anodic)',...
                'intBPE(i_cathodic)', ...
                '(intBPE(i_anodic)-intBPE(i_cathodic))/2', ... % faradaic current, averaged
                '(intWE(Ix)+intCE(Ix))/2 - (intBPE(i_anodic)-intBPE(i_cathodic))/2' }; % ohmic current, averaged
            
for i=1:numel(expressions)
    model.result.numerical.create(titles{i}, 'EvalGlobal');
    model.result.numerical(titles{i}).set('probetag', 'none');
    model.result.numerical(titles{i}).set('table', 'tbl1');
    model.result.numerical(titles{i}).set('expr', expressions{i});
    model.result.numerical(titles{i}).set('descr', expressions{i});
end

%% export
model.result.export.create('exportTertiaryCurrentDistributionData', 'Data');
% model.result.export('exportTertiaryCurrentDistributionData').set('data', dset);
model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', m.c_id{i}, i);
end

% files('exportTertiaryCurrentDistributionDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistributionDataFile.txt'];
model.result.export('exportTertiaryCurrentDistributionData').set('location', 'grid');
model.result.export('exportTertiaryCurrentDistributionData').set('gridx2', 'range( XleftBoundary, W/1000, XrightBoundary)');
model.result.export('exportTertiaryCurrentDistributionData').set('gridy2', '0');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', exportTertiaryCurrentDistributionDataFile);
% model.result.export('exportTertiaryCurrentDistributionData').run;

% full
model.result.export.create('exportTertiaryCurrentDistribution2dData', 'Data');
% model.result.export('exportTertiaryCurrentDistributionData').set('data', dset);
% model.result.export('exportTertiaryCurrentDistribution2dData').setIndex('expr', 'phi', 0);
% for i=1:m.numberOfSpecies
%     model.result.export('exportTertiaryCurrentDistribution2dData').setIndex('expr', m.c_id{i}, i);
% end
% 
% % files('exportTertiaryCurrentDistribution2dDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistribution2dDataFile.txt'];
% 
% model.result.export('exportTertiaryCurrentDistribution2dData').set('location', 'fromdataset');
% model.result.export('exportTertiaryCurrentDistribution2dData').set('filename', exportTertiaryCurrentDistribution2dDataFile);
% % model.result.export('exportTertiaryCurrentDistributionData').run;
% model.result.export('exportTertiaryCurrentDistributionData').set('location', 'fromdataset');
gridx = 10000;
% W = mphevaluate(model,'W');
% H = mphevaluate(model,'H');
% %     gridx = round(W/H*4);
% gridy = round(gridx*H/W);
gridy = round(gridx*m.H/m.W);

% model.result.export('exportTertiaryCurrentDistribution2dData').set('data', dset);
model.result.export('exportTertiaryCurrentDistribution2dData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistribution2dData').setIndex('expr', m.c_id{i}, i);
end

model.result.export('exportTertiaryCurrentDistribution2dData').set('location', 'regulargrid');
model.result.export('exportTertiaryCurrentDistribution2dData').set('resolution', 'finer');
model.result.export('exportTertiaryCurrentDistribution2dData').set('sort', 'on');
model.result.export('exportTertiaryCurrentDistribution2dData').set('regulargridx2', gridx);
model.result.export('exportTertiaryCurrentDistribution2dData').set('regulargridy2', gridy);
model.result.export('exportTertiaryCurrentDistribution2dData').set('filename', exportTertiaryCurrentDistribution2dDataFile);
model.result.export('exportTertiaryCurrentDistribution2dData').run;

m.saveState;