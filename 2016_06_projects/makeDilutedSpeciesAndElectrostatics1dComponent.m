% component for tertiary current distribution
%% create
% component for rough tertiary current approximation
% m.comp_id = 'tertiaryCurrentDistributionComponent';
model.modelNode.create('dilutedSpeciesAndElectrostatics1dComponent'); 
model.geom.create('dilutedSpeciesAndElectrostatics1dGeometry',1);
model.geom('dilutedSpeciesAndElectrostatics1dGeometry').insertFile(files('geometryPartsMphFile'), 'simple1dGeometry');

% mesh
model.mesh.create('dilutedSpeciesAndElectrostatics1dMesh', 'dilutedSpeciesAndElectrostatics1dGeometry');
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').create('copy1', 'Copy');

% since copy not possible for 1d, explicit:
m.hMaxFactor = 0.1;
m.mesh1D();

r = m.intFirstDebyeLength(2)/m.intFirstDebyeLength(1);
a0 = m.intFirstDebyeLength(1);

r_zeta = m.intFirstDebyeLength(end)/m.intFirstDebyeLength(end-1);
a_zeta = m.intExtendedDdl(1);

r_bulk = m.intRemaining(2)/m.intRemaining(1);
a0_bulk = m.intRemaining(1);

model.mesh('dilutedSpeciesAndElectrostatics1dMesh').create('edg1', 'Edge'); % ddl
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').create('dis1', 'Distribution');
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').create('size1', 'Size');

model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').selection.geom(1);
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').selection.all;

model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('dis1').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_ddl_dom');
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('dis1').set('type','explicit');
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('dis1').set('explicit', m.distributionFirstDebyeLengthStr);

model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_zetaVertex');
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('custom', 'on');
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hmaxactive', true);
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hgrad', r_zeta);
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hmax', sprintf('%e*L',a_zeta));
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hgradactive', true);

model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('size').set('custom', 'on');
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('size').set('hmin', sprintf('%e*L',a0));
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('size').set('hgrad', r_zeta);
model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('size').set('hmax', sprintf('%e*L',a0_bulk));

model.mesh('dilutedSpeciesAndElectrostatics1dMesh').run;

% functions
model.func.create('smoothenBpeBC', 'Rectangle');
model.func.create('interpolateStoredValues','Interpolation');

% analytical pb
% model.func.create('phi_pb', 'Analytic');
% model.func.create('phi_pbx', 'Analytic');
% for i = 1:m.numberOfSpecies   
%     model.func.create(m.c_pb_id{i}, 'Analytic');
% end
model.func.create('pbDecayFunction', 'Analytic');

%operators
model.cpl.create('onSurface', 'LinearExtrusion', 'dilutedSpeciesAndElectrostatics1dGeometry');
model.cpl.create('intSurface', 'Integration', 'dilutedSpeciesAndElectrostatics1dGeometry');
model.cpl.create('intBulk', 'Integration', 'dilutedSpeciesAndElectrostatics1dGeometry');


% variables
model.variable.create('domainVariables');
model.variable.create('surfaceVariables');
% model.variable.create('bulkVariables');

% physics
model.physics.create('DilutedSpecies', 'DilutedSpecies', 'dilutedSpeciesAndElectrostatics1dGeometry', m.c_id);
model.physics.create('Electrostatics', 'Electrostatics', 'dilutedSpeciesAndElectrostatics1dGeometry');

model.physics('DilutedSpecies').feature.create('BulkConcentration', 'Concentration', 0);
% model.physics('TertiaryCurrentDistribution').feature.create('BpeSurface', 'ExternalElectrodeSurface', 1);
model.physics('DilutedSpecies').feature.create('SurfaceFlux', 'Fluxes', 0);

model.physics('Electrostatics').feature.create('SpaceChargeDensity', 'SpaceChargeDensity', 1);
model.physics('Electrostatics').feature.create('SurfaceChargeDensity', 'SurfaceChargeDensity', 0);
model.physics('Electrostatics').feature.create('BulkPotential', 'ElectricPotential', 0);

% model.physics('TertiaryCurrentDistribution').feature.create('ElectrodePotential', 'ElectrolytePotential', 1);

% for i=1:m.nReactions
% %     model.physics('TertiaryCurrentDistribution').feature('BpeSurface').feature.create(m.reactionNames{i}, 'ElectrodeReaction', 1);
% end

%% update
% dummy parameter for sweep
model.param.set('X',0);

% dimensional space charge density term
chargeDensitySummand = prepTerm('z_id*c_id','z_id','c_id',m.z_id,m.c_id);
chargeDensityTerm = strcat('F_const*(',... % 
        strjoin(chargeDensitySummand','+' ), ')');
    
% dimensionless pb estimates
% phiPBFunction = '4/abs(z_ref)*atanh( tanh( abs(z_ref)*phi0/4 )*exp(-(x/epsilon+delta)))'; % dimensionless
% phiPBxFunction = '- 4/(epsilon*abs(z_ref))* sinh( abs(z_ref) * phi0/4 )* cosh( abs(z_ref) * phi0/4 ) / ( exp(x/epsilon+delta)* cosh( abs(z_ref) * phi0/4 )^2 - exp(-(x/epsilon+delta)) * sinh( abs(z_ref) * phi0/4 )^2)';

% dimensionless:
%   cPB* = c_bulk/c_ref * exp( -z*phiPB*)
% cPBFunction = prepTerm('c_bulk_id/c_ref*exp(-z_id*phi_pb(x,phi0))',...
%     'c_bulk_id','z_id',m.c_bulk_id,m.z_id);

pbDecayFunction = 'exp(-(x/lambdaD+delta))';

pdePSourceTerm = prepTerm('chargeDensityRampFactor*chargeDensityTerm','chargeDensityTerm',chargeDensityTerm); % allow for ramping


% mesh
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('copy1').set('mesh', 'simple1dGeometryRefinedMeshPart');
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('copy1').selection('source').geom(1);
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('copy1').selection('destination').geom(1);
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('copy1').selection('source').all;
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('copy1').selection('destination').all;
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').run('copy1');

% functions 
model.func('smoothenBpeBC').set('upper', 'w_bpe/2');
model.func('smoothenBpeBC').set('smooth', 'epsilon*smootheningFactor');
model.func('smoothenBpeBC').set('funcname', 'smoothenBpeBC');
model.func('smoothenBpeBC').set('lower', '-w_bpe/2');

% model.func('interpolateStoredValues').set('nargs', '2');
model.func('interpolateStoredValues').set('source', 'file');
model.func('interpolateStoredValues').set('interp', 'cubicspline');
model.func('interpolateStoredValues').setIndex('funcs', 'dummy_zero', 0, 0);
model.func('interpolateStoredValues').setIndex('funcs', num2str(1), 0, 1);

model.func('interpolateStoredValues').setIndex('funcs', 'phi_bulk_ddl', 1, 0);
model.func('interpolateStoredValues').setIndex('funcs', num2str(2), 1, 1);

for i=1:m.numberOfSpecies
    model.func('interpolateStoredValues').setIndex('funcs', sprintf('%s_bulk_ddl',m.c_id{i}), i+1, 0);
    model.func('interpolateStoredValues').setIndex('funcs', num2str(i+2), i+1, 1);
end
model.func('interpolateStoredValues').set('filename', files('exportTertiaryCurrentDistributionDataFile'));
model.func('interpolateStoredValues').set('nargs', '1');
model.func('interpolateStoredValues').importData;

% analytical pb
% model.func('phi_pb').set('funcname', 'phi_pb');
% model.func('phi_pb').set('args', {'x' 'phi0'});
% model.func('phi_pb').set('expr', phiPBFunction);
% model.func('phi_pb').set('fununit', '1');
% model.func('phi_pb').set('argunit', '1');
% model.func('phi_pb').set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
% model.func('phi_pb').label('Potential distribution in symmetric electrolyte due to Poisson-Boltzmann model');
% 
% model.func('phi_pbx').set('funcname', 'phi_pbx');
% model.func('phi_pbx').set('args', {'x' 'phi0'});
% model.func('phi_pbx').set('expr', phiPBxFunction);
% model.func('phi_pbx').set('fununit', '1');
% model.func('phi_pbx').set('argunit', '1');
% model.func('phi_pbx').set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
% model.func('phi_pbx').label('Derivative of potential distribution');
% 
% for i = 1:m.numberOfSpecies   
%     model.func(m.c_pb_id{i}).set('funcname', m.c_pb_id{i});
%     model.func(m.c_pb_id{i}).set('args', {'x' 'phi0'});
%     model.func(m.c_pb_id{i}).set('expr', cPBFunction{i});
%     model.func(m.c_pb_id{i}).set('fununit', '1');
%     model.func(m.c_pb_id{i}).set('argunit', '1');
%     model.func(m.c_pb_id{i}).set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
% end

model.func('pbDecayFunction').set('funcname', 'pbDecayFunction');
model.func('pbDecayFunction').set('args', {'x'});
model.func('pbDecayFunction').set('expr', pbDecayFunction);
model.func('pbDecayFunction').set('fununit', '1');
model.func('pbDecayFunction').set('argunit', 'm');
model.func('pbDecayFunction').set('plotargs', {'x' '0' 'L'});

% operators
model.cpl('onSurface').set('opname','onSurface');
% model.cpl('onSurface').selection('srcvertex2').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
% model.cpl('onSurface').selection('srcvertex1').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
% model.cpl('onSurface').selection('dstvertex2').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
% model.cpl('onSurface').selection('dstvertex1').geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
model.cpl('onSurface').selection.geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
model.cpl('onSurface').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
model.cpl('onSurface').selection('srcvertex1').set(1); % order of vertex creation important
model.cpl('onSurface').selection('dstvertex1').set(3);

model.cpl('intSurface').selection.geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
model.cpl('intSurface').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
model.cpl('intSurface').set('opname', 'intSurface');

model.cpl('intBulk').selection.geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
model.cpl('intBulk').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');
model.cpl('intBulk').set('opname', 'intBulk');

% variables
model.variable('domainVariables').model('dilutedSpeciesAndElectrostatics1dComponent');
model.variable('domainVariables').selection.geom('dilutedSpeciesAndElectrostatics1dGeometry', 1);
model.variable('domainVariables').selection.all;
model.variable('domainVariables').loadFile(files('domainVariablesDilutedSpeciesAndElectrostatics1d'));

model.variable('surfaceVariables').model('dilutedSpeciesAndElectrostatics1dComponent');
model.variable('surfaceVariables').selection.geom('dilutedSpeciesAndElectrostatics1dGeometry', 0);
model.variable('surfaceVariables').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
model.variable('surfaceVariables').loadFile(files('surfaceVariablesDilutedSpeciesAndElectrostatics1d'));


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


model.physics('Electrostatics').feature('SpaceChargeDensity').selection.all;
model.physics('Electrostatics').feature('SpaceChargeDensity').set('rhoq', pdePSourceTerm); % space charge density

% potential initial values
model.physics('Electrostatics').feature('init1').set('phi', 'phi0');

% bc selection
model.physics('DilutedSpecies').feature('SurfaceFlux').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
model.physics('DilutedSpecies').feature('BulkFlux').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');
model.physics('Electrostatics').feature('SurfaceChargeDensity').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');

model.physics('DilutedSpecies').feature('BulkConcentration').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');
model.physics('Electrostatics').feature('BulkPotential').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');

model.physics('Electrostatics').feature('SurfaceChargeDensity').set('rhoqs', 'smoothenBpeBC(X/L)*C_Stern_dimensional*(phi_s-phi_bulk_ddl(X))');
model.physics('Electrostatics').feature('BulkPotential').set('V0', 'phi_bulk_ddl(X)');


for i = 1:m.numberOfSpecies
    D_c_id = sprintf('D_%s',m.c_id{i});
    % isotropic diffusivity
    model.physics('DilutedSpecies').feature('cdm1').set(D_c_id, {m.D_id{i} '0' '0' '0' m.D_id{i} '0' '0' '0' m.D_id{i}});

    model.physics('DilutedSpecies').feature('cdm1').setIndex('z', m.z_id{i}, i-1);

    % bulk and initial concentrations
    
    model.physics('DilutedSpecies').feature('BulkConcentration').setIndex('species', true, i-1);
    model.physics('DilutedSpecies').feature('BulkConcentration').setIndex('c0', sprintf('%s_bulk_ddl(X)',m.c_id{i}), i-1);
    model.physics('DilutedSpecies').feature('init1').setIndex('initc', m.c_0_id{i}, i-1);

    
    model.physics('DilutedSpecies').feature('SurfaceFlux').setIndex('species', true, i-1);
    model.physics('DilutedSpecies').feature('SurfaceFlux').setIndex('N0', sprintf('smoothenBpeBC(X/L)*%s', m.N_id{i}), i-1);
 
    model.physics('DilutedSpecies').feature('BulkFlux').setIndex('species', true, i-1);
    model.physics('DilutedSpecies').feature('BulkFlux').setIndex('N0', sprintf('-smoothenBpeBC(X/L)*onSurface(%s)', m.N_id{i}), i-1);
 
end