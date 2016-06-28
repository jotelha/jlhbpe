% component for tertiary current distribution
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
model.physics.create('c', 'CoefficientFormPDE', 'dilutedSpeciesAndElectrostatics2dGeometry');

model.physics('c').feature.create('init2', 'init', 2);
% model.physics('Electrostatics').feature.create('init2', 'init', 2);

% model.physics('DilutedSpecies').feature.create('BulkConcentration', 'Concentration', 1);
model.physics('c').create('BulkDirichletBC', 'DirichletBoundary', 1);

% model.physics('DilutedSpecies').feature.create('DdlConcentration', 'Concentration', 1);
model.physics('c').create('DdlDirichletBC', 'DirichletBoundary', 1);

% model.physics('TertiaryCurrentDistribution').feature.create('BpeSurface', 'ExternalElectrodeSurface', 1);
% model.physics('DilutedSpecies').feature.create('SurfaceFlux', 'Fluxes', 1);
% model.physics('DilutedSpecies').feature.create('BulkFlux', 'Fluxes', 0);

% model.physics('Electrostatics').feature.create('SpaceChargeDensity', 'SpaceChargeDensity', 2);
% model.physics('Electrostatics').feature.create('SurfaceChargeDensity', 'SurfaceChargeDensity', 1);
% model.physics('Electrostatics').feature.create('BulkPotential', 'ElectricPotential', 1);
% model.physics('Electrostatics').feature.create('BulkPotential', 'ElectricPotential', 1);
% model.physics('Electrostatics').feature.create('DdlPotential', 'ElectricPotential', 1);

% model.physics('DilutedSpecies').create('DilutedSpeciesContinuity', 'Continuity', 1);
% model.physics('Electrostatics').create('ElectrostaticsContinuity', 'Continuity', 1);
model.physics('c').create('CoefficientPDEContinuity', 'Continuity', 1);

%% update

% dimensional space charge density term
chargeDensitySummand = prepTerm('z_id*c_id','z_id','c_id',m.z_id,m.c_id);
chargeDensityTerm = strcat('F_const*(',... % 
        strjoin(chargeDensitySummand','+' ), ')');
    
pbDecayFunction = 'exp(-(x/lambdaD+delta))';

pdePSourceTerm = prepTerm('chargeDensityRampFactor*chargeDensityTerm','chargeDensityTerm',chargeDensityTerm); % allow for ramping

pdeNPAlphaX = prepTerm('z_id*F_const/RT*D_id*phix','z_id','D_id',m.z_id,m.D_id);
pdeNPAlphaY = prepTerm('z_id*F_const/RT*D_id*phiy','z_id','D_id',m.z_id,m.D_id);


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
model.func('interpolateStoredValues1d').set('filename', files('exportDilutedSpeciesAndElectrostatics1dDataFile'));
model.func('interpolateStoredValues1d').set('nargs', '2');
model.func('interpolateStoredValues1d').importData;

% interpolate previous bulk results
model.func('interpolateStoredValues2d').set('source', 'file');
model.func('interpolateStoredValues2d').set('interp', 'linear');

model.func('interpolateStoredValues2d').setIndex('funcs', 'phi_interp_bulk', 0, 0);
model.func('interpolateStoredValues2d').setIndex('funcs', num2str(1), 0, 1);

for i=1:m.numberOfSpecies
    model.func('interpolateStoredValues2d').setIndex('funcs', sprintf('%s_interp_bulk',m.c_id{i}), i, 0);
    model.func('interpolateStoredValues2d').setIndex('funcs', num2str(i+1), i, 1);
end
model.func('interpolateStoredValues2d').set('filename', files('exportTertiaryCurrentDistribution2dDataFile'));
model.func('interpolateStoredValues2d').set('nargs', '2');
model.func('interpolateStoredValues2d').importData;


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

model.physics('c').field('dimensionless').field('c');
model.physics('c').field('dimensionless').component([m.c_id 'phi']);

model.physics('c').prop('ShapeProperty').set('order', '7');
model.physics('c').prop('ShapeProperty').set('valueType', 'real');


% standard settings
model.physics('c').feature('cfeq1').set('c', jlh.hf.flattenCell( jlh.hf.strdiag([m.D_id 'epsilon_r*epsilon0_const']) )' );
model.physics('c').feature('cfeq1').set('f', [jlh.hf.strzeros(m.numberOfSpecies,1); pdePSourceTerm]);
model.physics('c').feature('cfeq1').set('al',[ ...
    jlh.hf.flattenCell( jlh.hf.strdiag([pdeNPAlphaX; '0']) )', ...
    jlh.hf.flattenCell( jlh.hf.strdiag([pdeNPAlphaY; '0']) )' ] );

model.physics('c').feature('init2').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_simpleDdlGeometryPartInstance1_ddl_dom');

% potential initial values
model.physics('c').feature('init1').set('phi', 'phi_interp_bulk(x,y)');
model.physics('c').feature('init2').set('phi', 'phi_interp_ddl(lambdaD+y,x)');

% bc selection
model.physics('c').feature('BulkDirichletBC').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_bulkOuterBoundary');
model.physics('c').feature('DdlDirichletBC').selection.named('dilutedSpeciesAndElectrostatics2dGeometry_ddlOuterBoundary');

model.physics('c').feature('BulkDirichletBC').set('r', prepTerm('FuncName_interp_bulk(x,y)','FuncName',[m.c_id,'phi']) );
model.physics('c').feature('DdlDirichletBC').set('r', prepTerm('FuncName_interp_ddl(lambdaD+y,x)','FuncName',[m.c_id,'phi']) );

for i = 1:m.numberOfSpecies
    % bulk and initial concentrations
    model.physics('c').feature('init1').set(m.c_id{i}, sprintf('%s_interp_bulk(x,y)',m.c_id{i}));
    model.physics('c').feature('init2').set(m.c_id{i}, sprintf('%s_interp_ddl(lambdaD+y,x)',m.c_id{i}));
end

model.physics('c').feature('CoefficientPDEContinuity').setIndex('pairs', 'p1', 0);