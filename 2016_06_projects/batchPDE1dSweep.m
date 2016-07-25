import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end

caseStudyTitle = 'pde1d';
createEmptyProject;

% or, to load an existing project from server or from file
% tag = 'xyz'
% loadExistingProject

% prepare text files
% makeParameterFile
% makeVariablesFiles

if ~exist('hMaxFactor','var')
    hMaxFactor = 0.01;
end

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

% GlobalEvaluationsTableFileName = [OutputDir,'\globalEvaluations.txt'];
exportDilutedSpeciesAndElectrostatics1dDataFile = [OutputDir,'\exportDilutedSpeciesAndElectrostatics1dDataFile.txt'];
% exportTertiaryCurrentDistribution2dDataFile = [OutputDir,'\exportTertiaryCurrentDistribution2dDataFile.txt'];

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

% mesh parts

% model.modelNode.create('mcomp1', 'MeshComponent');
% 
% model.geom.create('mgeom1', 2); % mesh geometry
% model.mesh.create('simpleBulkGeometryRefinedMeshPart', 'mgeom1');
% model.mesh('simpleBulkGeometryRefinedMeshPart').create('imp1', 'Import');
% model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
% model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('filename', files('simpleBulkGeometryRefinedMeshFile'));
% model.mesh('simpleBulkGeometryRefinedMeshPart').run;
% model.mesh('mpart1').feature('imp1').set('facepartition', 'auto');

% model.modelNode.create('mcomp2', 'MeshComponent');
% model.geom.create('mgeom2', 2); % mesh geometry
% model.mesh.create('simpleDdlGeometryRefinedMeshPart', 'mgeom2');
% model.mesh('simpleDdlGeometryRefinedMeshPart').create('imp1', 'Import');
% model.mesh('simpleDdlGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
% model.mesh('simpleDdlGeometryRefinedMeshPart').feature('imp1').set('filename', files('simpleDdlGeometryRefinedMeshFile'));
% model.mesh('simpleDdlGeometryRefinedMeshPart').run;

%% makeTertiaryCurrentDistribution2dComponentWithNetCurrentConst

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
m.hMaxFactor = hMaxFactor;
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

% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_zetaVertex');
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('custom', 'on');
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hmaxactive', true);
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hgrad', r_zeta);
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hmax', sprintf('%e*L',a_zeta));
% model.mesh('dilutedSpeciesAndElectrostatics1dMesh').feature('edg1').feature('size1').set('hgradactive', true);

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
% model.physics.create('DilutedSpecies', 'DilutedSpecies', 'dilutedSpeciesAndElectrostatics1dGeometry', m.c_id);
% model.physics.create('Electrostatics', 'Electrostatics', 'dilutedSpeciesAndElectrostatics1dGeometry');
model.physics.create('c', 'CoefficientFormPDE', 'dilutedSpeciesAndElectrostatics1dGeometry');

% model.physics('DilutedSpecies').feature.create('BulkConcentration', 'Concentration', 0);
% model.physics('TertiaryCurrentDistribution').feature.create('BpeSurface', 'ExternalElectrodeSurface', 1);
% model.physics('DilutedSpecies').feature.create('SurfaceFlux', 'Fluxes', 0);
% model.physics('DilutedSpecies').feature.create('BulkFlux', 'Fluxes', 0);
model.physics('c').create('SurfaceFlux', 'FluxBoundary', 0);
model.physics('c').create('BulkDirichlet', 'DirichletBoundary', 0);


% model.physics('Electrostatics').feature.create('SpaceChargeDensity', 'SpaceChargeDensity', 1);
% model.physics('Electrostatics').feature.create('SurfaceChargeDensity', 'SurfaceChargeDensity', 0);
% model.physics('Electrostatics').feature.create('BulkPotential', 'ElectricPotential', 0);

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

pdeNPAlpha = prepTerm('z_id*F_const/RT*D_id*phix','z_id','D_id',m.z_id,m.D_id);

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
% model.func('interpolateStoredValues').set('filename', files('exportTertiaryCurrentDistributionDataFile'));
% model.func('interpolateStoredValues').set('nargs', '1');
% model.func('interpolateStoredValues').importData;

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
% model.cpl('onSurface').selection('dstvertex1').set(3);
model.cpl('onSurface').selection('dstvertex1').set(2);

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

% model.physics('DilutedSpecies').prop('TransportMechanism').set('Convection', false);
% model.physics('DilutedSpecies').prop('TransportMechanism').set('Migration', true);
% 
% model.physics('Electrostatics').field('electricpotential').field('phi');
% model.physics('DilutedSpecies').field('concentration').component(m.c_id);
model.physics('c').field('dimensionless').field('c');
% c1 to cN on 1st to Nth row, phi on N+1th row
model.physics('c').field('dimensionless').component([m.c_id 'phi']);

% model.physics('Electrostatics').prop('ShapeProperty').set('order_electricpotential', '5');
% model.physics('Electrostatics').prop('ShapeProperty').set('valueType', 'real');
model.physics('c').prop('ShapeProperty').set('order', '7');
model.physics('c').prop('ShapeProperty').set('valueType', 'real');

% model.physics('DilutedSpecies').prop('ShapeProperty').set('order_concentration', '5');
% model.physics('DilutedSpecies').prop('ShapeProperty').set('valueType', 'real');

% standard settings
% model.physics('DilutedSpecies').feature('cdm1').set('minput_temperature', 'T');
% model.physics('DilutedSpecies').feature('cdm1').set('V', 'phi');
% model.physics('Electrostatics').feature('ccn1').set('epsilonr_mat', 'userdef');
% model.physics('Electrostatics').feature('ccn1').set('epsilonr', {'epsilon_r' '0' '0' '0' 'epsilon_r' '0' '0' '0' 'epsilon_r'});

model.physics('c').feature('cfeq1').set('c', jlh.hf.flattenCell( jlh.hf.strdiag([m.D_id 'epsilon_r*epsilon0_const']) )' );
model.physics('c').feature('cfeq1').set('f', [jlh.hf.strzeros(m.numberOfSpecies,1); pdePSourceTerm]);
model.physics('c').feature('cfeq1').set('al',jlh.hf.flattenCell( jlh.hf.strdiag([pdeNPAlpha; '0']) ) );

% model.physics('Electrostatics').feature('SpaceChargeDensity').selection.all;
% model.physics('Electrostatics').feature('SpaceChargeDensity').set('rhoq', pdePSourceTerm); % space charge density

% potential initial values
% model.physics('Electrostatics').feature('init1').set('phi', 'phi0');
model.physics('c').feature('init1').set('phi', 'phi0');

% bc selection
% model.physics('DilutedSpecies').feature('SurfaceFlux').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
% model.physics('DilutedSpecies').feature('BulkFlux').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');
% model.physics('Electrostatics').feature('SurfaceChargeDensity').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
model.physics('c').feature('SurfaceFlux').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_surfaceVertex');
model.physics('c').feature('SurfaceFlux').set('g', prepTerm('smoothenBpeBC(X/L)*FluxBC','FluxBC',[m.N_id 'C_Stern_dimensional*phi_s']));
model.physics('c').feature('SurfaceFlux').set('q', jlh.hf.flattenCell( jlh.hf.strdiag( [jlh.hf.strzeros(1,m.numberOfSpecies) 'smoothenBpeBC(X/L)*C_Stern_dimensional'] ) ) );

% model.physics('DilutedSpecies').feature('BulkConcentration').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');
% model.physics('Electrostatics').feature('BulkPotential').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');
model.physics('c').feature('BulkDirichlet').selection.named('dilutedSpeciesAndElectrostatics1dGeometry_simple1dGeometryPartInstance1_bulkVertex');

% model.physics('Electrostatics').feature('SurfaceChargeDensity').set('rhoqs', 'smoothenBpeBC(X/L)*C_Stern_dimensional*(phi_s-phi_bulk_ddl(X))');
% model.physics('Electrostatics').feature('BulkPotential').set('V0', 'phi_bulk_ddl(X)');
model.physics('c').feature('BulkDirichlet').set('r', prepTerm('FuncName_bulk_ddl(X)','FuncName',[m.c_id,'phi']) );

for i = 1:m.numberOfSpecies
    % bulk and initial concentrations
    
%     model.physics('DilutedSpecies').feature('BulkConcentration').setIndex('species', true, i-1);
%     model.physics('DilutedSpecies').feature('BulkConcentration').setIndex('c0', sprintf('%s_bulk_ddl(X)',m.c_id{i}), i-1);
%     model.physics('DilutedSpecies').feature('init1').setIndex('initc', m.c_0_id{i}, i-1);
    model.physics('c').feature('init1').set(m.c_id{i}, m.c_0_id{i});

    
%     model.physics('DilutedSpecies').feature('SurfaceFlux').setIndex('species', true, i-1);
%     model.physics('DilutedSpecies').feature('SurfaceFlux').setIndex('N0', sprintf('smoothenBpeBC(X/L)*%s', m.N_id{i}), i-1);
%  
%     model.physics('DilutedSpecies').feature('BulkFlux').setIndex('species', true, i-1);
%     model.physics('DilutedSpecies').feature('BulkFlux').setIndex('N0', sprintf('-smoothenBpeBC(X/L)*onSurface(%s)', m.N_id{i}), i-1);
%  
end

%% meshing
% if( ~exist('meshFile','var' ) )
%     meshFile = 'copySimpleBulkGeometryMeshCoarse_case_a.m'; %% standard method of meshing, copy mesh from geometry file
% end
% 
% run(meshFile);

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
model.batch('p1').feature('cl1').set('filename', 'UpdateInitialValues1d.class');
model.batch('p1').feature('cl1').set('input', {'sweep' 'input\inputFiles.txt' 'input' 'interpolateStoredValues'});
model.batch(param_id).feature('so1').set('seq', sol_id);
model.batch('p1').feature('nu1').set('table', 'tbl1');

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

% for x-coordinate sweep
inner_pname = {'X'};
inner_plistarr = {'range(XleftBoundary,W/1000,XrightBoundary)'};

model.study(study_id).feature(studyStep_id).set('useparam', 'on');
model.study(study_id).feature(studyStep_id).set('plistarr', inner_plistarr);
model.study(study_id).feature(studyStep_id).set('pname', inner_pname);
model.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
model.study(study_id).feature(studyStep_id).set('pcontinuationmode', 'no');
model.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
model.study(study_id).feature(studyStep_id).set('showdistribute', true);

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

% for x-coordinate sweep
model.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
model.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', 'no');
model.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', inner_plistarr);
model.sol(sol_id).feature(solver_id).feature('pDef').set('pname', inner_pname); 

% 2d parametric extrusion
% make parametric extrusion
model.result.dataset.create('par1', 'Parametric1D');
% model.result.dataset('par1').set('data', dset);
% model.result.dataset('par1').setIndex('looplevelinput', 'manual', 2);
% model.result.dataset('par1').setIndex('looplevelinput', 'first', 1);
% model.result.dataset('par1').setIndex('looplevelinput', 'manual', 1);
model.result.dataset('par1').set('levelscaleactive', 'on');
model.result.dataset('par1').set('levelscale', '1');
model.result.dataset('par1').set('innerinput', 'all');
% model.result.dataset('par1').run;

% %% numerical evaluations
% % model.result.table('tbl1').comments('Global Evaluation 1 (intWE(tcdee.Ilx))');
% model.result.export.create('tbl1', 'Table');
% model.result.export('tbl1').set('filename', GlobalEvaluationsTableFileName);
% 
% titles = {'DeltaPHI','PHI_bpe', 'Ix_WE', 'Ix_CE', 'Iy_BPE', ' I_total', 'I_anodic', 'I_cathodic', 'I_faradaic', 'I_ohmic'};
% expressions = { 'DeltaPHI', ...
%                 'intWE(tcdee.Ilx)', ...
%                 'intCE(tcdee.Ilx)', ...
%                 'intBPE(tcdee.Ily)', ...
%                 'intBPE(i_total)',...
%                 'intBPE(i_anodic)',...
%                 'intBPE(i_cathodic)', ...
%                 '(intBPE(i_anodic)-intBPE(i_cathodic))/2', ... % faradaic current, averaged
%                 '(intWE(tcdee.Ilx)+intCE(tcdee.Ilx))/2 - (intBPE(i_anodic)-intBPE(i_cathodic))/2' }; % ohmic current, averaged
%             
% for i=1:numel(expressions)
%     model.result.numerical.create(titles{i}, 'EvalGlobal');
%     model.result.numerical(titles{i}).set('probetag', 'none');
%     model.result.numerical(titles{i}).set('table', 'tbl1');
%     model.result.numerical(titles{i}).set('expr', expressions{i});
%     model.result.numerical(titles{i}).set('descr', expressions{i});
% end

%% export
model.result.export.create('exportTertiaryCurrentDistributionData', 'Data');
model.result.export('exportTertiaryCurrentDistributionData').set('data', 'par1');
model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', m.c_id{i}, i);
end

% files('exportTertiaryCurrentDistributionDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistributionDataFile.txt'];
model.result.export('exportTertiaryCurrentDistributionData').set('location', 'fromdataset');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', exportDilutedSpeciesAndElectrostatics1dDataFile);
% model.result.export('exportTertiaryCurrentDistributionData').run;


m.saveState;