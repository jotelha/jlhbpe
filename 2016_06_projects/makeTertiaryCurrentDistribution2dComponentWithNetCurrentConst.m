% component for tertiary current distribution
%% create
% component for rough tertiary current approximation
% m.comp_id = 'tertiaryCurrentDistributionComponent';
model.modelNode.create('tertiaryCurrentDistributionComponent'); 
model.geom.create('tertiaryCurrentDistributionGeometry',2);
model.geom('tertiaryCurrentDistributionGeometry').insertFile(files('geometryPartsMphFile'), 'simpleBulkGeometry');

%model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
%model.mesh('tertiaryCurrentDistributionMesh').create('copy1', 'Copy');

% model.geom('tertiaryCurrentDistributionGeometry').insertFile(geometryPartsMphFile, 'simpleAssembledGeometry');

% functions
model.func.create('smoothenBpeBC', 'Rectangle');
model.func.create('interpolateStoredValues2d','Interpolation');

% operators
model.cpl.create('intWE', 'Integration', 'tertiaryCurrentDistributionGeometry');
model.cpl.create('intCE', 'Integration', 'tertiaryCurrentDistributionGeometry');
model.cpl.create('intBPE', 'Integration', 'tertiaryCurrentDistributionGeometry');

% variables
model.variable.create('domainVariables');
model.variable.create('surfaceVariables');
model.variable.create('weVariables');
model.variable.create('ceVariables');

% physics
model.physics.create('TertiaryCurrentDistribution', 'TertiaryCurrentDistributionNernstPlanck', 'tertiaryCurrentDistributionGeometry', m.c_id);
model.physics('TertiaryCurrentDistribution').feature.create('BulkConcentration', 'Concentration', 1);
model.physics('TertiaryCurrentDistribution').feature.create('BpeSurface', 'ExternalElectrodeSurface', 1);
model.physics('TertiaryCurrentDistribution').feature.create('ElectrodePotential', 'ElectrolytePotential', 1);

for i=1:m.nReactions
    model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature.create(m.reactionNames{i}, 'ElectrodeReaction', 1);
end

model.physics.create('ge', 'GlobalEquations', 'tertiaryCurrentDistributionGeometry');
model.physics('ge').feature.create('gconstr1', 'GlobalConstraint', -1);
%% update

% mesh
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').geom(2);
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').geom(2);
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').all;
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_space_dom');
% model.mesh('tertiaryCurrentDistributionMesh').run('copy1');
run(meshFile);

% functions 
model.func('smoothenBpeBC').set('upper', 'w_bpe/2');
model.func('smoothenBpeBC').set('smooth', 'epsilon*smootheningFactor');
model.func('smoothenBpeBC').set('funcname', 'smoothenBpeBC');
model.func('smoothenBpeBC').set('lower', '-w_bpe/2');

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

% operators
model.cpl('intWE').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_workingElectrode');
model.cpl('intCE').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_counterElectrode');
model.cpl('intBPE').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_bpeSurface');

model.cpl('intWE').label('intWE');
model.cpl('intWE').set('opname', 'intWE');
model.cpl('intCE').label('intCE');
model.cpl('intCE').set('opname', 'intCE');
model.cpl('intBPE').label('intBPE');
model.cpl('intBPE').set('opname', 'intBPE');

% variables
model.variable('domainVariables').model('tertiaryCurrentDistributionComponent');
model.variable('domainVariables').selection.geom('tertiaryCurrentDistributionGeometry', 2);
model.variable('domainVariables').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_space_dom');
model.variable('domainVariables').loadFile(files('domainVariablesTertiaryCurrentDistribution2d'));

model.variable('surfaceVariables').model('tertiaryCurrentDistributionComponent');
model.variable('surfaceVariables').selection.geom('tertiaryCurrentDistributionGeometry', 1);
model.variable('surfaceVariables').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_bpeSurface');
model.variable('surfaceVariables').loadFile(files('surfaceVariablesTertiaryCurrentDistribution2d'));

model.variable('weVariables').model('tertiaryCurrentDistributionComponent');
model.variable('weVariables').selection.geom('tertiaryCurrentDistributionGeometry', 1);
model.variable('weVariables').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_workingElectrode');
model.variable('weVariables').loadFile(files('weVariablesTertiaryCurrentDistribution2d'));

model.variable('ceVariables').model('tertiaryCurrentDistributionComponent');
model.variable('ceVariables').selection.geom('tertiaryCurrentDistributionGeometry', 1);
model.variable('ceVariables').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_counterElectrode');
model.variable('ceVariables').loadFile(files('ceVariablesTertiaryCurrentDistribution2d'));

% physics
model.physics('TertiaryCurrentDistribution').field('electricpotentialionicphase').field('phi');
model.physics('TertiaryCurrentDistribution').field('electricpotential').field('phi_e'); % redundant
model.physics('TertiaryCurrentDistribution').field('concentration').component(m.c_id);
model.physics('TertiaryCurrentDistribution').prop('ShapeProperty').set('order_concentration', '5');
model.physics('TertiaryCurrentDistribution').prop('ShapeProperty').set('order_electricpotentialionicphase', '5');
model.physics('TertiaryCurrentDistribution').prop('ShapeProperty').setIndex('valueType', 'real', 0, 0);
model.physics('TertiaryCurrentDistribution').prop('ShapeProperty').setIndex('valueType', 'real', 2, 0);
model.physics('TertiaryCurrentDistribution').prop('ShapeProperty').setIndex('valueType', 'real', 1, 0);

% last species from electroneutrality
model.physics('TertiaryCurrentDistribution').prop('SpeciesProperties').set('FromElectroneutrality', m.numberOfSpecies);

% ice1 is the standard electrolyte settings
model.physics('TertiaryCurrentDistribution').feature('ice1').set('minput_temperature', 'T');

% potential initial values
model.physics('TertiaryCurrentDistribution').feature('init1').set('initphil', 'phi_interp_bulk(x,y)');
model.physics('TertiaryCurrentDistribution').feature('init1').set('initphis', 'PHI_bpe'); % redundant

% bc selection
model.physics('TertiaryCurrentDistribution').feature('BpeSurface').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_bpeSurface');
model.physics('TertiaryCurrentDistribution').feature('ElectrodePotential').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_electrodes');
% model.physics('TertiaryCurrentDistribution').feature('BulkConcentration').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_bulkBoundary');
model.physics('TertiaryCurrentDistribution').feature('BulkConcentration').selection.named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_electrodes');

model.physics('TertiaryCurrentDistribution').feature('BpeSurface').set('phisext0', 'PHI_bpe'); % bpe
model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature('er1').active(false); % deactivate standard reaction
model.physics('TertiaryCurrentDistribution').feature('ElectrodePotential').set('philbnd', 'phi_s'); % feeder electrodes

for i = 1:m.numberOfSpecies
    D_c_id = sprintf('D_%s',m.c_id{i});
    % isotropic diffusivity
    model.physics('TertiaryCurrentDistribution').feature('ice1').set(D_c_id, {m.D_id{i} '0' '0' '0' m.D_id{i} '0' '0' '0' m.D_id{i}});

    model.physics('TertiaryCurrentDistribution').feature('ice1').setIndex('z', m.z_id{i}, i-1);

    % bulk and initial concentrations
    if i<m.numberOfSpecies
        model.physics('TertiaryCurrentDistribution').feature('BulkConcentration').setIndex('species', true, i-1);
        model.physics('TertiaryCurrentDistribution').feature('BulkConcentration').setIndex('c0', m.c_bulk_id{i}, i-1);
        model.physics('TertiaryCurrentDistribution').feature('init1').setIndex('initc', sprintf('%s_interp_bulk(x,y)',m.c_id{i}), i-1);
    end
end
    
for i=1:m.nReactions
    model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature(m.reactionNames{i}).setIndex('minput_temperature', 'T', 0);
    model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature(m.reactionNames{i}).set('ElectrodeKinetics', 'userdef');
    model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature(m.reactionNames{i}).set('iloc', sprintf('smoothenBpeBC(x/L)*%s',m.i_id{i}));
    model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature(m.reactionNames{i}).set('nm', m.n_id{i});
    for j =1:(m.numberOfSpecies-1) % stochiometric coefficients
        model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature(m.reactionNames{i}).setIndex('Vi0', m.nu_id{j,i}, j-1);
    end
end

% zero net current on bpe surface equation:
model.physics('ge').feature('ge1').set('name', 'PHI_bpe');
model.physics('ge').feature('ge1').set('equation', 'intBPE(tcdee.Ily)');
model.physics('ge').feature('ge1').set('initialValueU', m.PHI_bpe);
model.physics('ge').feature('ge1').set('valueType', 'real');

% constraint for stabilization: entering current must equal exiting current
model.physics('ge').feature('gconstr1').set('constraintExpression', 'intWE(tcdee.Ily)-intCE(tcdee.Ily)');