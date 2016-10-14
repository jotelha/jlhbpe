import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end

caseStudyTitle = 'tcd2d';
createEmptyProject;

% or, to load an existing project from server or from file
% tag = 'xyz'
% loadExistingProject

% prepare text files
% makeParameterFile
% makeVariablesFiles

if ~exist('BatchDir','var')
    BatchDir = 'batch';
end
if ~exist('ParametricDir','var')
    ParametricDir = 'parametric';
end
if ~exist('OutputDir','var')
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

%% makeTertiaryCurrentDistribution2dComponent
% component for rough tertiary current approximation
model.modelNode.create('tertiaryCurrentDistributionComponent'); 
model.geom.create('tertiaryCurrentDistributionGeometry',2);
model.geom('tertiaryCurrentDistributionGeometry').insertFile(files('geometryPartsMphFile'), 'simpleBulkGeometry');

% model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
% model.mesh('tertiaryCurrentDistributionMesh').create('copy1', 'Copy');

% model.geom('tertiaryCurrentDistributionGeometry').insertFile(geometryPartsMphFile, 'simpleAssembledGeometry');

% functions
model.func.create('smoothenBpeBC', 'Rectangle');

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

%% update

% mesh
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').geom(2);
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').geom(2);
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').all;
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_space_dom');
% model.mesh('tertiaryCurrentDistributionMesh').run('copy1');

% functions 
model.func('smoothenBpeBC').set('upper', 'w_bpe/2');
model.func('smoothenBpeBC').set('smooth', 'epsilon*smootheningFactor');
model.func('smoothenBpeBC').set('funcname', 'smoothenBpeBC');
model.func('smoothenBpeBC').set('lower', '-w_bpe/2');

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
model.physics('TertiaryCurrentDistribution').feature('init1').set('initphil', 'phi0');
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
        model.physics('TertiaryCurrentDistribution').feature('init1').setIndex('initc', m.c_0_id{i}, i-1);
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

% makeBatchStudy
% tags of created features in sol_id, study_id, studyStep_id, compile_id,
% variables_id, solver_id and store_id
% model.sol('sol1').feature('s1').feature('dDef').set('ooc', 'on');

% %% meshing
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

model.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
% model.study(study_id).feature(studyStep_id).set('showdistribute', true);

% model.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);
model.sol(sol_id).feature(solver_id).set('probesel', 'none');
model.sol(sol_id).feature(solver_id).set('nonlin', 'on');
model.sol(sol_id).feature(solver_id).set('stol', '1e-3');

% model.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'const');
model.sol(sol_id).feature(solver_id).feature('dDef').set('ooc', 'off'); % out of core
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
                'PHI_bpe', ... % mixed potential
                'intWE(tcdee.Ilx)', ...
                'intCE(tcdee.Ilx)', ...
                'intBPE(tcdee.Ily)', ...
                'intBPE(i_total)',...
                'intBPE(i_anodic)',...
                'intBPE(i_cathodic)', ...
                '(intBPE(i_anodic)-intBPE(i_cathodic))/2', ... % faradaic current, averaged
                '(intWE(tcdee.Ilx)+intCE(tcdee.Ilx))/2 - (intBPE(i_anodic)-intBPE(i_cathodic))/2' }; % ohmic current, averaged
            
for i=1:numel(expressions)
    model.result.numerical.create(titles{i}, 'EvalGlobal');
    model.result.numerical(titles{i}).set('probetag', 'none');
    model.result.numerical(titles{i}).set('table', 'tbl1');
    model.result.numerical(titles{i}).set('expr', expressions{i});
    model.result.numerical(titles{i}).set('descr', expressions{i});
end



%% export
%% export
model.result.export.create('exportTertiaryCurrentDistributionData', 'Data');
model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', m.c_id{i}, i);
end

model.result.export('exportTertiaryCurrentDistributionData').set('location', 'grid');
model.result.export('exportTertiaryCurrentDistributionData').set('gridx2', 'range( XleftBoundary, W/1000, XrightBoundary)');
model.result.export('exportTertiaryCurrentDistributionData').set('gridy2', '0');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', exportTertiaryCurrentDistributionDataFile);

% full
model.result.export.create('exportTertiaryCurrentDistribution2dData', 'Data');

gridx = 1000;
gridy = round(gridx*m.H/m.W);

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