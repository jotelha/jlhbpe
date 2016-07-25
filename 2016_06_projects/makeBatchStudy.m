%% single-step bpe study, 2d
% sweepParameterValues = {-0.28210/m.UT};
% sweepParameterValues = {-0.01/m.UT};

if ~exist('BatchDir','var')
%     BatchDir = [ pwd, '\', m.projectPath, '\batch'];
    BatchDir = 'batch';
end
if ~exist('ParametricDir','var')
    ParametricDir = 'parametric';
end
BatchDirAbsolute = [pwd(), '\', m.projectPath, '\' BatchDir];
ParametricDirAbsolute = [pwd(), '\', m.projectPath, '\' ParametricDir];

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
model.batch('p1').feature('cl1').set('filename', 'UpdateInitialValues.class');
model.batch('p1').feature('cl1').set('input', {'sweep' '"testList.txt"' '".\input"' '"int1"'});
model.batch(param_id).feature('so1').set('seq', sol_id);
model.batch('p1').feature('nu1').set('table', 'tbl1');

model.batch(param_id).attach(study_id);

% batch
model.study(study_id).create(batchStep_id, 'Batch');
model.batch.create(batch_id,'Batch');

model.batch(batch_id).create('jo1', 'Jobseq');
% model.batch(batch_id).feature('daDef').create('pr1', 'Process');

% model.batch(batch_id).create('nu1', 'Numericalseq');
% model.batch(batch_id).create('ex1', 'Exportseq');
model.batch(batch_id).study(study_id);

model.study(study_id).feature(batchStep_id).set('batchdir', BatchDir);
model.study(study_id).feature(batchStep_id).set('batchfile', BatchFileName);

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

% what do those stand for:
% model.study(study_id).feature(studyStep_id).set('initstudyhide', 'on');
% model.study(study_id).feature(studyStep_id).set('initsolhide', 'on');
% model.study(study_id).feature(studyStep_id).set('solnumhide', 'on');
% model.study(study_id).feature(studyStep_id).set('notstudyhide', 'on');
% model.study(study_id).feature(studyStep_id).set('notsolhide', 'on');
% model.study(study_id).feature(studyStep_id).set('notsolnumhide', 'on');

% model.study(study_id).feature(studyStep_id).set('mesh', {'geom' 'mappedMesh' 'geom1d' 'nomesh'});

% model.study(study_id).feature(studyStep_id).set('useparam', useparam);
% model.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
% model.study(study_id).feature(studyStep_id).set('pname', pname);
% model.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
% model.study(study_id).feature(studyStep_id).set('pcontinuationmode', pcontinuationmode);
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
% model.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
% model.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
% model.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
% model.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
% model.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname); 
