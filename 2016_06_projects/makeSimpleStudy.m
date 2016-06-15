%% single-step bpe study, 2d
% sweepParameterValues = {-0.28210/m.UT};
% sweepParameterValues = {-0.01/m.UT};

% pname = {'phi_bpe_init','phi_bpe','surfaceFluxRampFactor', 'deltaPhi'};
% plistarr = cellstr( cellfun(@(c) num2str(c), sweepParameterValues,'UniformOutput',false));
pname = {'deltaPhi','surfaceFluxRampFactor'};
% pval = cellstr( cellfun(@(c) num2str(c), sweepPhiInitValues,'UniformOutput',false));
plistarr = ['1'; '1'];

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

model.study(study_id).create(studyStep_id, 'Stationary');

model.sol(sol_id).create(compile_id, 'StudyStep');
model.sol(sol_id).create(variables_id, 'Variables');
model.sol(sol_id).create(solver_id, 'Stationary');
model.sol(sol_id).create(store_id, 'StoreSolution');

% set scaling for variables
model.sol(sol_id).feature(variables_id).set('scalemethod','init');

% what do those stand for:
model.study(study_id).feature(studyStep_id).set('initstudyhide', 'on');
model.study(study_id).feature(studyStep_id).set('initsolhide', 'on');
model.study(study_id).feature(studyStep_id).set('solnumhide', 'on');
model.study(study_id).feature(studyStep_id).set('notstudyhide', 'on');
model.study(study_id).feature(studyStep_id).set('notsolhide', 'on');
model.study(study_id).feature(studyStep_id).set('notsolnumhide', 'on');

% model.study(study_id).feature(studyStep_id).set('mesh', {'geom' 'mappedMesh' 'geom1d' 'nomesh'});

model.study(study_id).feature(studyStep_id).set('useparam', useparam);
model.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
model.study(study_id).feature(studyStep_id).set('pname', pname);
model.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
model.study(study_id).feature(studyStep_id).set('pcontinuationmode', pcontinuationmode);
model.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
model.study(study_id).feature(studyStep_id).set('showdistribute', true);

model.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);
model.sol(sol_id).feature(solver_id).set('probesel', 'none');
model.sol(sol_id).feature(solver_id).set('nonlin', 'on');
model.sol(sol_id).feature(solver_id).set('stol', '1e-3');

% model.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'const');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
% model.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermconst', 'itertol');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '15');
% model.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
% model.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
model.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
model.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
model.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
model.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname); 
