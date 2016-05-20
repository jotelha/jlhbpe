 import com.comsol.model.*
 import com.comsol.model.util.*
 import jlh.*
    
 %%
m = jlh.BpeModel;
m.newProject('Duval2003FaradaicSmall1d');
%% also necessary to execute when loading
% server log file
logFile = [pwd(),'\',m.projectPath,'\comsol.log'];
% fclose(fopen(logFile, 'w'));
ModelUtil.showProgress(logFile);
system(sprintf('G:\\scripts\\mtail\\mTail.exe "%s" /start &',logFile));

% m.setDhopeshwarkar2008ElectrokineticsCaseParameters;

setDuval2003FaradaicParametersSmall

% E_m = m.calcMixedPotential;
% m.plotCurrents;
m.PHI_bpe = E_m(c);

m.prepareIdentifiers;

%% m.createModel;
loadedModels    = ModelUtil.tags;
isLoaded        = arrayfun( @(s) strcmp(m.model_tag,s),loadedModels);
if any(isLoaded)
    ModelUtil.remove(m.model_tag);
end
m.m = ModelUtil.create(m.model_tag);   % creates model on COMSOL server
m.m.disableUpdates(true); % try whether speed can be increased

m.createFunctions()
m.m.modelNode.create(m.comp_id); 
m.updateParameters();
m.makeChoppedGeometry();

%% 
% m.createGeometry();
% m.createSelections();
m.createOperators()
m.createVariables();
m.createPhysicsForWeakForm();
m.createMesh();
m.createProbes();
m.create1dPlots();
m.create2dPlots();
m.standardView = m.m.view.create('standardView', 'geom');
m.ddlView      = m.m.view.create('ddlView', 'geom');
m.createDatasets();

m.create0dComponent();

%% 1d component
m.create1dComponent();
m.createEvaluations1d;
m.update1dGeometry();
m.update1dOperators();
%% m.updateModel;

m.updateFunctions();
% m.updateSelections();
% m.updateOperators();
m.updateOperatorsForChoppedGeometry();
% m.updatePhysicsForChoppedGeometry();
m.updatePhysicsForWeakForm();

% m.updateMeshMapped();
% m.addParametricSweepStudy('chargeDensityRampFactor',0);
% m.addParametricSweepStudy({'chargeDensityRampFactor','phiSymmetryFactor'},{[0,1]; [0:0.1:1]});
% m.addParametricSweepStudy({'deltaPhiRampFactor'},{1});
% m.setUniformFeederElectrodePotentialBC
% m.setUniformFeederElectrodePotentialAndLinearBulkPotentialBC
% m.setFeederElectrodeCurrentFluxAndBPEPotentialBC
% m.setMinimalBC;
% m.setUniformFeederElectrodePotentialAndLinearBulkPotentialBC;
% m.activateBpe;
m.m.disableUpdates(false); % necessary before meshing

m.saveState;

%% meshing for 1d component
m.hMaxFactor = 0.01;
m.mesh1D;
m.update1dMesh();

m.saveState;
%% study for 1d component
sweepParameterValues = {-0.28210*m.UT};
pname = {'phi_bpe_init'};
plistarr = cellstr( cellfun(@(c) num2str(c), sweepParameterValues,'UniformOutput',false));
punit = '';
useparam = 'on';
pcontinuationmode = 'no';

          
studyNum = m.m.study.size() + 1;
study_id = sprintf('study%u',studyNum);

currentStationaryStudy = m.m.study.create(study_id);

solNum = m.m.sol.size() + 1;
sol_id = sprintf('sol%u',solNum);

m.m.sol.create(sol_id);
m.m.sol(sol_id).study(study_id); % ?
m.m.sol(sol_id).attach(study_id); % ?


% for i=1:nStudySteps
i = 1;
fprintf('\t Creating step %d and setting standard parameters...\n',i);
studyStep_id = sprintf('stationaryStudyStep%d',i);
compile_id = sprintf('st%d',i);
variables_id = sprintf('v%d',i);
solver_id = sprintf('s%d',i);
store_id = sprintf('su%d',i);

m.m.study(study_id).create(studyStep_id, 'Stationary');

m.m.sol(sol_id).create(compile_id, 'StudyStep');
m.m.sol(sol_id).create(variables_id, 'Variables');
m.m.sol(sol_id).create(solver_id, 'Stationary');
m.m.sol(sol_id).create(store_id, 'StoreSolution');

% what do those stand for:
m.m.study(study_id).feature(studyStep_id).set('initstudyhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('initsolhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('solnumhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('notstudyhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('notsolhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('notsolnumhide', 'on');

m.m.study(study_id).feature(studyStep_id).set('mesh', {'geom' 'nomesh' 'geom1d' 'explicitMesh1d'});

m.m.study(study_id).feature(studyStep_id).set('useparam', useparam);
m.m.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
m.m.study(study_id).feature(studyStep_id).set('pname', pname);
m.m.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
m.m.study(study_id).feature(studyStep_id).set('pcontinuationmode', pcontinuationmode);
m.m.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
m.m.study(study_id).feature(studyStep_id).set('showdistribute', true);


m.m.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);
% 
%     % standard consecutive solving order:
%     if( i > 1)
%         m.m.study(study_id).feature(studyStep_id).set('solnum', 'auto');
%         m.m.study(study_id).feature(studyStep_id).set('useinitsol', 'on');
%         m.m.study(study_id).feature(studyStep_id).set('initmethod', 'sol');
%         m.m.study(study_id).feature(studyStep_id).set('initstudy', study_id);
% 
%         m.m.sol(sol_id).feature(variables_id).set('initsol', sol_id);
%         m.m.sol(sol_id).feature(variables_id).set('solnum', 'auto');
%         m.m.sol(sol_id).feature(variables_id).set('initmethod', 'sol');
%     end

m.m.sol(sol_id).feature(solver_id).set('probesel', 'none');
m.m.sol(sol_id).feature(solver_id).set('nonlin', 'on');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '30');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname); 

display('Setting individual active physics for step...');

m.m.study(study_id).feature(studyStep_id).set('disabledphysics', {''});
m.m.study(study_id).feature(studyStep_id).activate('WeakFormulation', false);
% m.m.study(study_id).feature(studyStep_id).activate('WeakFormulation1d', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroSurfaceCurrent', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroNetCurrent0d', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroNetCurrent1d', false);

m.saveState;
%% solve 1d 
m.m.sol(sol_id).run(store_id);

%% plot 1d plots
storedSolutions = m.m.study(study_id).getSolverSequences('Stored');
solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(1)));
dset = char(solInfo.dataset(2));

m.updateEvaluations1d(dset);
m.plot1dSolution(dset);

%% meshing for 2d component

% 
% m.hMaxFactor = 0.1;
% m.updateChoppedMesh();
% m.replicateMeshPrototype;
% m.finalizeChoppedMesh;
% m.standardMesh.run;
% m.saveState;

% mapped mesh
m.m.disableUpdates(false); % necessary before meshing
m.hMaxFactor = 1;
m.createMappedMesh();
m.updateMappedMesh();
m.replicateMappedMeshPrototype;
% m.finalizeMappedMesh;
m.m.mesh('mappedMesh').run;
m.saveState;
m.m.disableUpdates(true); 

% m.mesh1D;
% a = m.intFirstDebyeLength(1);
% r = m.intFirstDebyeLength(2)/m.intFirstDebyeLength(1);
% N = ceil(log(1-m.epsilon/a*(1-r))/log(r));
%% multi-step bpe study
% nSteps = 5;
nStudySteps = 5;
mesh_id = 'mappedMesh';
sweepParameterValues = {1};

fprintf('Creating multi-step bpe study with %d steps...\n',nStudySteps);
% m.addBpeStudyStepSequence({'deltaPhiRampFactor'},{1},...
%     {   {   'PoissonEquation', ...
%             'NernstPlanckEquation',...
%             'zeroSurfaceCurrent'} ,...
%         {   'NernstPlanckEquation/FluxAtSurface', ...
%             'zeroSurfaceCurrent',...
%             'zeroNetCurrent0d'}, ...
%         {   'PoissonEquation', ...
%             'NernstPlanckEquation',...
%             'zeroNetCurrent0d' }, ...
%          {  'zeroSurfaceCurrent',...
%             'zeroNetCurrent0d'}  ,...
%          {  'zeroNetCurrent0d'}  });
%         


pname = {'deltaPhiRampFactor'};
plistarr = cellstr( cellfun(@(c) num2str(c), sweepParameterValues,'UniformOutput',false));
punit = '';
useparam = 'on';
pcontinuationmode = 'no';
    % OR
%     pcontinuationmode = 'manual'; % or 'auto', 'last'
%     pcontinuation = pname;
%     pout = 'plist'; % or 'psteps'

          
studyNum = m.m.study.size() + 1;
study_id = sprintf('study%u',studyNum);

stationaryStudy = m.m.study.create(study_id);

solNum = m.m.sol.size() + 1;
sol_id = sprintf('sol%u',solNum);

m.m.sol.create(sol_id);
m.m.sol(sol_id).study(study_id); % ?
m.m.sol(sol_id).attach(study_id); % ?


for i=1:nStudySteps
    fprintf('\t Creating step %d and setting standard parameters...\n',i);
    studyStep_id = sprintf('stationaryStudyStep%d',i);
    compile_id = sprintf('st%d',i);
    variables_id = sprintf('v%d',i);
    solver_id = sprintf('s%d',i);
    store_id = sprintf('su%d',i);

    stationaryStudy.create(studyStep_id, 'Stationary');

    m.m.sol(sol_id).create(compile_id, 'StudyStep');
    m.m.sol(sol_id).create(variables_id, 'Variables');
    m.m.sol(sol_id).create(solver_id, 'Stationary');
    m.m.sol(sol_id).create(store_id, 'StoreSolution');

%         m.m.sol(sol_id).feature(solver_id).feature.remove('fcDef');
%         m.m.sol(sol_id).feature(solver_id).create('fc1', 'FullyCoupled');
%         m.m.sol(sol_id).feature(solver_id).create('pDef', 'Parametric');
%         m.m.sol(sol_id).feature(solver_id).create('p1', 'Parametric');

    % what do those stand for:
    m.m.study(study_id).feature(studyStep_id).set('initstudyhide', 'on');
    m.m.study(study_id).feature(studyStep_id).set('initsolhide', 'on');
    m.m.study(study_id).feature(studyStep_id).set('solnumhide', 'on');
    m.m.study(study_id).feature(studyStep_id).set('notstudyhide', 'on');
    m.m.study(study_id).feature(studyStep_id).set('notsolhide', 'on');
    m.m.study(study_id).feature(studyStep_id).set('notsolnumhide', 'on');

%     m.m.study(study_id).feature(studyStep_id).set('disabledphysics', disablePhysics{i});
    m.m.study(study_id).feature(studyStep_id).set('mesh', {'geom' mesh_id 'geom1d' 'nomesh'});

    m.m.study(study_id).feature(studyStep_id).set('useparam', useparam);
    m.m.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
    m.m.study(study_id).feature(studyStep_id).set('pname', pname);
    m.m.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
    m.m.study(study_id).feature(studyStep_id).set('pcontinuationmode', pcontinuationmode);
    m.m.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
    m.m.study(study_id).feature(studyStep_id).set('showdistribute', true);


    m.m.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);

    % standard consecutive solving order:
    if( i > 1)
        m.m.study(study_id).feature(studyStep_id).set('solnum', 'auto');
        m.m.study(study_id).feature(studyStep_id).set('useinitsol', 'on');
        m.m.study(study_id).feature(studyStep_id).set('initmethod', 'sol');
        m.m.study(study_id).feature(studyStep_id).set('initstudy', study_id);

        m.m.sol(sol_id).feature(variables_id).set('initsol', sol_id);
        m.m.sol(sol_id).feature(variables_id).set('solnum', 'auto');
        m.m.sol(sol_id).feature(variables_id).set('initmethod', 'sol');
    end


    m.m.sol(sol_id).feature(solver_id).set('probesel', 'none');
    m.m.sol(sol_id).feature(solver_id).set('nonlin', 'on');
    m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
    m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
    m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '30');
    m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
    m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
    m.m.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
    m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
    m.m.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
    m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname);
end 

%% set the 'correct' disable mode
display('Setting individual active physics for each step...');

% m.m.study('study1').feature('stationaryStudyStep1').set('disabledphysics', {'NernstPlanckEquation' 'PoissonEquation' 'zeroSurfaceCurrent'});
% m.m.study('study1').feature('stationaryStudyStep1').set('disabledphysics', {''});
m.m.study('study1').feature('stationaryStudyStep1').set('disabledphysics', {''});
% m.m.study('study1').feature('stationaryStudyStep1').activate('NernstPlanckEquation', false);
m.m.study('study1').feature('stationaryStudyStep1').activate('WeakFormulation', false);
% m.m.study('study1').feature('stationaryStudyStep1').activate('PoissonEquation', false);
m.m.study('study1').feature('stationaryStudyStep1').activate('WeakFormulation1d', false);
m.m.study('study1').feature('stationaryStudyStep1').activate('zeroSurfaceCurrent', false);
m.m.study('study1').feature('stationaryStudyStep1').activate('zeroNetCurrent1d', false);

% m.m.study('study1').feature('stationaryStudyStep2').set('disabledphysics', {'NernstPlanckEquation/FluxAtSurface' 'zeroSurfaceCurrent' 'zeroNetCurrent0d'});
% m.m.study('study1').feature('stationaryStudyStep2').set('disabledphysics', {'NernstPlanckEquation/FluxAtSurface'});
m.m.study('study1').feature('stationaryStudyStep2').set('pname', [pname,{'surfaceFluxRampFactor'}]);
m.m.study('study1').feature('stationaryStudyStep2').set('plistarr', [plistarr,'0']);
m.m.study('study1').feature('stationaryStudyStep2').activate('WeakFormulation1d', false);
m.m.study('study1').feature('stationaryStudyStep2').activate('zeroSurfaceCurrent', false);
m.m.study('study1').feature('stationaryStudyStep2').activate('zeroNetCurrent0d', false);
m.m.study('study1').feature('stationaryStudyStep2').activate('zeroNetCurrent1d', false);

% m.m.study('study1').feature('stationaryStudyStep3').set('disabledphysics', {'NernstPlanckEquation' 'PoissonEquation' 'zeroNetCurrent0d'});
m.m.study('study1').feature('stationaryStudyStep3').set('disabledphysics', {''});
% m.m.study('study1').feature('stationaryStudyStep3').activate('NernstPlanckEquation', false);
% m.m.study('study1').feature('stationaryStudyStep3').activate('PoissonEquation', false);
m.m.study('study1').feature('stationaryStudyStep3').activate('WeakFormulation', false);
m.m.study('study1').feature('stationaryStudyStep3').activate('WeakFormulation1d', false);
m.m.study('study1').feature('stationaryStudyStep3').activate('zeroSurfaceCurrent', false);
m.m.study('study1').feature('stationaryStudyStep3').activate('zeroNetCurrent1d', false);


% m.m.study('study1').feature('stationaryStudyStep4').set('disabledphysics', {'zeroSurfaceCurrent' 'zeroNetCurrent0d'});
m.m.study('study1').feature('stationaryStudyStep4').set('disabledphysics', {''});
m.m.study('study1').feature('stationaryStudyStep4').activate('WeakFormulation1d', false);
m.m.study('study1').feature('stationaryStudyStep4').activate('zeroSurfaceCurrent', false);
m.m.study('study1').feature('stationaryStudyStep4').activate('zeroNetCurrent0d', false);
m.m.study('study1').feature('stationaryStudyStep4').activate('zeroNetCurrent1d', false);

% m.m.study('study1').feature('stationaryStudyStep4').set('disabledphysics', {'zeroSurfaceCurrent' 'zeroNetCurrent0d'});
m.m.study('study1').feature('stationaryStudyStep5').set('disabledphysics', {''});
m.m.study('study1').feature('stationaryStudyStep5').activate('zeroNetCurrent0d', false);
m.m.study('study1').feature('stationaryStudyStep5').activate('WeakFormulation1d', false);
m.m.study('study1').feature('stationaryStudyStep5').activate('zeroNetCurrent1d', false);


m.m.study('study1').setStoreSolution(true);

storedSolutions = m.m.study('study1').getSolverSequences('Stored');

m.saveState;

display('Setting individual initial and fixed values for each step...');

%% study step 1
%nothing

%% study step 2

% set initial value source
% m.m.study('study1').feature('stationaryStudyStep2').set('solnum', 'auto');
m.m.study('study1').feature('stationaryStudyStep2').set('initmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep2').set('useinitsol', 'on');
m.m.study('study1').feature('stationaryStudyStep2').set('initmethod', 'init');
m.m.study('study1').feature('stationaryStudyStep2').set('initstudy', 'zero');

% set values for variables not solved for
m.m.study('study1').feature('stationaryStudyStep2').set('usesol', 'on');
m.m.study('study1').feature('stationaryStudyStep2').set('notsolmethod', 'init');
m.m.study('study1').feature('stationaryStudyStep2').set('notstudy', 'zero');

% m.m.sol('sol1').run('su2');

% solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(2)));
% dset = char(solInfo.dataset);
% m.updateDatasets(dset);
% m.update1dPlots;
% m.updateSurfacePlots(dset);
% m.update2dPlots(dset);


%% study step 3

% set initial value source
m.m.study('study1').feature('stationaryStudyStep3').set('initmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep3').set('initstudy', 'study1');
m.m.study('study1').feature('stationaryStudyStep3').set('initsol', 'sol1');
m.m.study('study1').feature('stationaryStudyStep3').set('initsoluse', char(storedSolutions(2)));

% set values for unsolved
m.m.study('study1').feature('stationaryStudyStep3').set('usesol', 'on');
m.m.study('study1').feature('stationaryStudyStep3').set('notsolmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep3').set('notstudy', 'study1');
m.m.study('study1').feature('stationaryStudyStep3').set('notsol', 'sol1');
m.m.study('study1').feature('stationaryStudyStep3').set('notsoluse', char(storedSolutions(2)));


% m.PHI_bpe = mphglobal(m.m,'comp0d.V','dataset',dset);
% m.m.param.set('phi_bpe_init',m.phi_bpe,'dimensionless potential at surface of bpe');

%% study step 4

% set initial value source
m.m.study('study1').feature('stationaryStudyStep4').set('initmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep4').set('initstudy', 'study1');
m.m.study('study1').feature('stationaryStudyStep4').set('initsol', 'sol1');
m.m.study('study1').feature('stationaryStudyStep4').set('initsoluse', char(storedSolutions(3)));

% set values for unsolved
m.m.study('study1').feature('stationaryStudyStep4').set('usesol', 'on');
m.m.study('study1').feature('stationaryStudyStep4').set('notsolmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep4').set('notstudy', 'study1');
m.m.study('study1').feature('stationaryStudyStep4').set('notsol', 'sol1');
m.m.study('study1').feature('stationaryStudyStep4').set('notsoluse', char(storedSolutions(3)));

% m.m.sol('sol1').runFromTo('st3','su4');
% solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(4)));
% dset = char(solInfo.dataset);


%% study step 5

% set initial value source
m.m.study('study1').feature('stationaryStudyStep5').set('initmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep5').set('initstudy', 'study1');
m.m.study('study1').feature('stationaryStudyStep5').set('initsol', 'sol1');
m.m.study('study1').feature('stationaryStudyStep5').set('initsoluse', char(storedSolutions(4)));

% set values for unsolved
m.m.study('study1').feature('stationaryStudyStep5').set('usesol', 'on');
m.m.study('study1').feature('stationaryStudyStep5').set('notsolmethod', 'sol');
m.m.study('study1').feature('stationaryStudyStep5').set('notstudy', 'study1');
m.m.study('study1').feature('stationaryStudyStep5').set('notsol', 'sol1');
m.m.study('study1').feature('stationaryStudyStep5').set('notsoluse', char(storedSolutions(4)));

m.saveState;

%% run all study steps
m.m.disableUpdates(false);  % necessary before solving

display('\tSolving step 1...');
m.m.sol('sol1').run('su1');

solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(1)));
dset = char(solInfo.dataset);
fprintf('\tModifying BPE potential intial value from %f.3V ', m.PHI_bpe);
m.PHI_bpe = mphglobal(m.m,'comp0d.V','dataset',dset);
fprintf('to %f.3V...\n', m.PHI_bpe);
m.m.param.set('phi_bpe_init',m.phi_bpe,'dimensionless potential at surface of bpe');

display('\tSolving step 2...');
m.m.sol('sol1').runFromTo('st2','su2');
display('\tSolving step 3...');
m.m.sol('sol1').runFromTo('st3','su3');

solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(3)));
dset = char(solInfo.dataset);
phi_bpe = mphglobal(m.m,'comp1.phi_bpe','dataset',dset);
fprintf('\tModifying BPE potential intial value from %f.3V ', m.PHI_bpe);
m.PHI_bpe = phi_bpe*m.UT;
fprintf('to %f.3V...\n', m.PHI_bpe);
m.m.param.set('phi_bpe_init',m.phi_bpe,'dimensionless potential at surface of bpe');

display('\tSolving step 4...');
m.m.sol('sol1').runFromTo('st4','su4');

display('\tSolving step 5...');
m.m.sol('sol1').runFromTo('st5','su5');
solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(5)));
dset = char(solInfo.dataset);
phi_bpe = mphglobal(m.m,'comp1.phi_bpe','dataset',dset);
m.PHI_bpe = phi_bpe*m.UT;
fprintf('\tBPE now at %f.3V.\n', m.PHI_bpe );

m.saveState;

%% output plots
nSteps = 5;
storedSolutions = m.m.study('study1').getSolverSequences('Stored');
delay = 3;
for i = 1:nSteps
    solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(i)));
    dset = char(solInfo.dataset);
    while true
        try
            fprintf('Trying to save surface flux plots for %s...\n',dset');
            m.updateSurfacePlots(dset);
            break;
        catch
            fprintf('Failed, retrying in %d seconds...\n',delay);
            pause(delay);
        end
    end
        
    switch i
        case {1,3} ;
        case {2,4}
%             while true
%                 try
%                     fprintf('Trying to update datasets for %s...\n',dset');
%                     m.updateDatasets(dset);
%                     break;
%                 catch
%                     fprintf('Failed, retrying in %d seconds...\n',delay);
%                     pause(delay);
%                 end
%             end
            while true
                try
                    fprintf('Trying to save 1d plots for %s...\n',dset');                    
%                     m.update1dPlots;
                    m.plotPresetCrossections(dset);
                    break;
                catch
                    fprintf('Failed, retrying in %d seconds...\n',delay);                    
                    pause(delay);
                end
            end
            while true
                try
                    fprintf('Trying to save 2d plots for %s...\n',dset');                    
                    m.update2dPlots(dset);
                    break;
                catch
                    fprintf('Failed, retrying in %d seconds...\n',delay);                    
                    pause(delay);
                end
            end
                
          case {5}
            while true
                try
                    fprintf('Trying to save 1d preset crossection plots for %s...\n',dset');                    
%                     m.update1dPlots;
                    m.plotPresetCrossections(dset);
                    break;
                catch
                    fprintf('Failed, retrying in %d seconds...\n',delay);                    
                    pause(delay);
                end
            end
            while true
                try
                    fprintf('Trying to save 1d horizontal crossection sweep plots for %s...\n',dset');                    
                    m.sweepHorizontalCrossection(dset,0.05)
                    break;
                catch
                    fprintf('Failed, retrying in %d seconds...\n',delay);                    
                    pause(delay);
                end
            end
            while true
                try
                    fprintf('Trying to save 1d vertical crossection sweep plots for %s...\n',dset');                    
                    m.sweepVerticalCrossection(dset,m.w/10,0.3);
                    break;
                catch
                    fprintf('Failed, retrying in %d seconds...\n',delay);                    
                    pause(delay);
                end
            end
            while true
                try
                    fprintf('Trying to save 2d plots for %s...\n',dset');                    
                    m.update2dPlots(dset);
                    break;
                catch
                    fprintf('Failed, retrying in %d seconds...\n',delay);                    
                    pause(delay);
                end
            end
    end
    
end
    
% m.update2dPlots(dset)
% %% run
% m.stationaryStudy.run
% 
% %% a second study?
% m.addBpeStudy({'deltaPhiRampFactor'},{1});
% phi_bpe = -30:1:70;
% m.addBpeStudy({'phi_bpe'},{phi_bpe});

%% plot last sol / dset

storedSolutions = m.m.study('study1').getSolverSequences('Stored');
solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(end)));
dset = char(solInfo.dataset);
delay = 3;
while true
    try
        fprintf('Trying to save 1d preset crossection plots for %s...\n',dset');                    
        m.plotPresetCrossections(dset);
        break;
    catch
        fprintf('Failed, retrying in %d seconds...\n',delay);                    
        pause(delay);
    end
end
while true
    try
        fprintf('Trying to save 1d horizontal crossection sweep plots for %s...\n',dset');                    
        m.sweepHorizontalCrossection(dset,0.05)
        break;
    catch
        fprintf('Failed, retrying in %d seconds...\n',delay);                    
        pause(delay);
    end
end
while true
    try
        fprintf('Trying to save 1d vertical crossection sweep plots for %s...\n',dset');                    
        m.sweepVerticalCrossection(dset,m.w/10,0.3);
        break;
    catch
        fprintf('Failed, retrying in %d seconds...\n',delay);                    
        pause(delay);
    end
end
while true
    try
        fprintf('Trying to save 2d plots for %s...\n',dset');                    
        m.update2dPlots(dset);
        break;
    catch
        fprintf('Failed, retrying in %d seconds...\n',delay);                    
        pause(delay);
    end
end

%% magnification
while true
    try
        fprintf('Trying to save 1d horizontal crossection sweep plots for %s...\n',dset');                    
        m.sweepHorizontalCrossection(dset,0.005,'-w_bpe/2-2*epsilon','-w_bpe/2+2*epsilon',0,0.03)
        break;
    catch
        fprintf('Failed, retrying in %d seconds...\n',delay);                    
        pause(delay);
    end
end
%% after getting initial values or result
% % dset = m.getLatestDataset();
% 
% m.updateDatasets(dset)
% % m.iterateStandardPlots
% m.update1dPlots
% m.updateSurfacePlots(dset)
% % save 2d plots
% % m.updateViews
% m.update2dPlots(dset)
% 
% %% parametric plots
% m.updateGlobalPlotsBySolnum(dset,'phiSymmetryFactor');
% m.update1dParametricPlots(dset);
% m.update2dParametricPlots(dset);
% 

return

%% load from server
m = jlh.BpeModel;

tag = '2016_05_20_13_02_40_Duval2003FaradaicSmall1d';

m.projectName = tag;
m.projectPath = ['dat\',tag];
m.model_tag = tag;

% mphFile = [m.projectPath,'\',tag,'.mph'];
reloadModelFile(m)


setDuval2003FaradaicParametersSmall
% E_m = m.calcMixedPotential;
% m.plotCurrents;
m.PHI_bpe = E_m(c);
m.prepareIdentifiers;

% tags = ModelUtil.tags;
% m.loadFromServer( char( tags(1) ));
% mphload(mphFile,tag);
% m.loadFromServer( char( tag ));
