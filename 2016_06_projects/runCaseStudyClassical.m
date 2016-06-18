import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

% parameter file for geometry to be created
caseStudyParameterFile = 'parameters_duval2001bipolar.m';


caseStudyTitle = 'duval2001bipolar';
createEmptyProject;

% or, to load an existing project from server or from file
% tag = 'xyz'
% loadExistingProject

% prepare text files
% makeParameterFile
% makeVariablesFiles

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

%%
makeTertiaryCurrentDistribution2dComponent

makeSimpleStudy
% tags of created features in sol_id, study_id, studyStep_id, compile_id,
% variables_id, solver_id and store_id

model.sol(sol_id).run(store_id);

storedSolutions = model.study(study_id).getSolverSequences('Stored');
solInfo = mphsolinfo(model,'soltag',char(storedSolutions(1)));
dset = char(solInfo.dataset);

createPlots

%%
% set up expressions to plot
% small letters: dimensionles, capital letters: dimensional

logc = prepTerm('log(c_id)','c_id',m.c_id);
% cx = prepTerm('c_idx','c_id',obj.c_id);
% cy = prepTerm('c_idy','c_id',obj.c_id);
cx = prepTerm('tcdee.grad_c_idx','c_id',m.c_id);
cy = prepTerm('tcdee.grad_c_idy','c_id',m.c_id);
absc  = prepTerm('sqrt(cx^2+cy^2)','cx','cy',cx,cy);
% dimensional
C = prepTerm('c_id/c_ref','c_id',m.c_id);
logC = prepTerm('log(c_id/c_ref)','c_id',m.c_id);
Cx = prepTerm('(tcdee.grad_c_idx/c_ref*L)','c_id',m.c_id);
Cy = prepTerm('(tcdee.grad_c_idy/c_ref*L)','c_id',m.c_id);
absC  = prepTerm('sqrt(Cx^2+Cy^2)','Cx','Cy',Cx,Cy);

absn = prepTerm('sqrt(nx_id^2+ny_id^2)','nx_id','ny_id',m.nx_id,m.ny_id);
absN = prepTerm('sqrt(Nx_id^2+Ny_id^2)','Nx_id','Ny_id',m.Nx_id,m.Ny_id);

absDiffusiveFlux = prepTerm('tcdee.dfluxMag_c_id','c_id',m.c_id);
diffusiveFluxX = prepTerm('tcdee.dflux_c_idx','c_id',m.c_id);
diffusiveFluxY = prepTerm('tcdee.dflux_c_idy','c_id',m.c_id);
electrophoreticFluxX = prepTerm('tcdee.mflux_c_idx','c_id',m.c_id);
electrophoreticFluxY = prepTerm('tcdee.mflux_c_idy','c_id',m.c_id);
absElectrophoreticFlux = prepTerm('sqrt(efx^2+efy^2)','efx','efy',electrophoreticFluxX,electrophoreticFluxY);
absTotalFlux = prepTerm('tcdee.tfluxMag_c_id','c_id',m.c_id);
totalFluxX = prepTerm('tcdee.tflux_c_idx','c_id',m.c_id);
totalFluxY = prepTerm('tcdee.tflux_c_idy','c_id',m.c_id);
       
            
plots = containers.Map;
% plots(title) = { {expression1, expression2, ...}, ylabel };
plots('phi')    = { {'phi/UT'}, 'phi / U_T'};
plots('phix')   = { {'phix/UT*L'}, 'phi_x * L / U_T'};
plots('phiy')   = { {'phiy/UT*L'}, 'phi_y * L / U_T'};
plots('_PHI')   = { {'phi'}, 'phi / V' };
plots('_PHIx')  = { {'phix'}, 'phi_x / V m^{-1}'};
plots('_PHIy')  = { {'phiy'}, 'phi_y / V m^{-1}'};

plots('c')      = { C, 'c / c_ref' };
plots('logc')   = { logC, 'log(c/c_ref)'};
plots('absc')   = { absC,'|grad(c)| * L / c_ref'}; 
plots('cx')     = { Cx,'c_x * L / c_ref'};
plots('cy')     = { Cy, 'c_y * L / c_ref'};
plots('_C')     = { m.c_id, 'c / mol m^{-3}'};
plots('_logC')  = { logc , 'log(c/mol m^{-3})'};
plots('_absC')  = { absc, '|grad c| / mol m^{-4}'};
plots('_Cx')    = { cx,  'c_x / mol m^{-4}'};
plots('_Cy')    = { cy, 'c_y / mol m^{-4}'};

plots('absn')   = { absn, '|n|'};
plots('nx')     = { m.nx_id,'nx'};
plots('ny')     = { m.ny_id,'ny'};
plots('absN')   = { absN, '|N|'};
plots('Nx')     = { m.Nx_id,'Nx'};
plots('Ny')     = { m.Ny_id,'Ny'};

plots('absDiffusiveFlux') = { absDiffusiveFlux, '-D*|grad c|' };
plots('diffusiveFluxX') = { diffusiveFluxX, '-D*cx' };
plots('diffusiveFluxY') = { diffusiveFluxY, '-D*cy' };

plots('absElectrophoreticFlux') = { absElectrophoreticFlux, '|-z*u*c*F*grad phi|' };
plots('electrophoreticFluxX') = { electrophoreticFluxX, '-z*u*c*F*phix' };
plots('electrophoreticFluxY') = { electrophoreticFluxY, '-z*u*c*F*phiy' };

plots('absTotalFlux') = { absTotalFlux, '|-D*grad c - z*u*c*F*grad phi|' };
plots('totalFluxX') = { totalFluxX, '-D*cx - z*u*c*F*phix' };
plots('totalFluxY') = { totalFluxY, '-z*u*c*F*phiy' };


plots('absi')   = { {'sqrt(ix^2+iy^2)'},'|i|'};
plots('ix')     = { {'ix'},'ix'};
plots('iy')     = { {'iy'},'iy'};
% plots('absI')   = { {'sqrt(Ix^2+iy^2)'},'|I|'};
% plots('Ix')     = { {'Ix'},'Ix'};
% plots('Iy')     = { {'Iy'},'Iy'};
plots('absI')   = { {'tcdee.IlMag'},'|I|'};
plots('Ix')     = { {'tcdee.Ilx'},'Ix'};
plots('Iy')     = { {'tcdee.Ily'},'Iy'};

plots('kappa')  = { {'kappa'},'kappa'};

sweepHorizontalCrossection(m,dset,m.L/4,'XleftBoundary','XrightBoundary',0,m.L,plots);
% m.sweepVerticalCrossection(dset,m.W/4);

%% export
model.result.export.create('exportTertiaryCurrentDistributionData', 'Data');
model.result.export('exportTertiaryCurrentDistributionData').set('data', dset);
model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', m.c_id{i}, i);
end

files('exportTertiaryCurrentDistributionDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistributionDataFile.txt'];

model.result.export('exportTertiaryCurrentDistributionData').set('location', 'grid');
model.result.export('exportTertiaryCurrentDistributionData').set('gridx2', 'range( XleftBoundary, W/1000, XrightBoundary)');
model.result.export('exportTertiaryCurrentDistributionData').set('gridy2', '0');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', files('exportTertiaryCurrentDistributionDataFile'));
model.result.export('exportTertiaryCurrentDistributionData').run;

return

%% geometry sequences

% % component for rough tertiary current approximation
% % m.comp_id = 'tertiaryCurrentDistributionComponent';
% model.modelNode.create('tertiaryCurrentDistributionComponent'); 
% model.geom.create('tertiaryCurrentDistributionGeometry',2);
% % model.geom('tertiaryCurrentDistributionGeometry').insertFile(geometryPartsMphFile, 'simpleAssembledGeometry');
% model.geom('tertiaryCurrentDistributionGeometry').insertFile(geometryPartsMphFile, 'simpleBulkGeometry');
% 
% model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
% model.mesh('tertiaryCurrentDistributionMesh').create('copy1', 'Copy');
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').geom(2);
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').geom(2);
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').all;
% model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_space_dom');
% model.mesh('tertiaryCurrentDistributionMesh').run('copy1');

% component for exact PNP solution
model.modelNode.create('dilutedSpeciesAndElectrostaticsComponent'); 
model.geom.create('dilutedSpeciesAndElectrostaticsGeometry',2);
% model.geom('tertiaryCurrentDistributionGeometry').insertFile(geometryPartsMphFile, 'simpleAssembledGeometry');
model.geom('dilutedSpeciesAndElectrostaticsGeometry').insertFile(files('geometryPartsMphFile'), 'simpleAssembledGeometry');

model.mesh.create('dilutedSpeciesAndElectrostaticsMesh', 'dilutedSpeciesAndElectrostaticsGeometry');

model.mesh('dilutedSpeciesAndElectrostaticsMesh').create('copy1', 'Copy');
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy1').selection('source').geom(2);
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy1').selection('destination').geom(2);
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy1').selection('source').all;
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy1').selection('destination').named('dilutedSpeciesAndElectrostaticsGeometry_simpleBulkGeometryPartInstance1_space_dom');
model.mesh('dilutedSpeciesAndElectrostaticsMesh').run('copy1');

model.mesh('dilutedSpeciesAndElectrostaticsMesh').create('copy2', 'Copy');
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy2').set('mesh', 'simpleDdlGeometryRefinedMeshPart');
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy2').selection('source').geom(2);
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy2').selection('destination').geom(2);
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy2').selection('source').all;
model.mesh('dilutedSpeciesAndElectrostaticsMesh').feature('copy2').selection('destination').named('dilutedSpeciesAndElectrostaticsGeometry_simpleDdlGeometryPartInstance1_ddl_dom');
model.mesh('dilutedSpeciesAndElectrostaticsMesh').run('copy2');

return;

% model.mesh(
% model.geom('tertiaryCurrentDistributionGeometry').create('tertiaryCurrentDistributionGeometryImport','Import');
% model.geom('tertiaryCurrentDistributionGeometry').feature('tertiaryCurrentDistributionGeometryImport').
%% create other features

m.createFunctions()
m.updateParameters();


m.createOperators()
m.createVariables();
% m.createPhysicsForWeakForm();
m.createPhysicsForClassicalForm();
% m.createMesh();
% m.createProbes();
m.create1dPlots();
m.create2dPlots();

m.m.result.export.create('dataExporter', 'Data');

m.standardView = m.m.view.create('standardView', 'geom');
m.ddlView      = m.m.view.create('ddlView', 'geom');
m.createDatasets();

m.create0dComponent();

%% 1d component
m.create1dComponent();
m.createEvaluations1d;
m.update1dGeometry();
m.update1dOperators();

%% set model characteristics
m.updateFunctions();
m.updateOperatorsForChoppedGeometry();
m.updatePhysicsForClassicalForm();

% m.m.disableUpdates(false); % necessary before meshing

m.saveState;

%% meshing for 1d component
m.hMaxFactor = 0.01;
m.mesh1D;
m.update1dMesh();

m.saveState;
%% study for 1d component
% sweepPhiInitValues = {-0.28210/m.UT};

% sweepPhiInitValues = {-0.01/m.UT};
% 
% pname = {'phi_bpe_init','surfaceFluxRampFactor'};
% plistarr = [ cellstr( cellfun(@(c) num2str(c), sweepPhiInitValues,'UniformOutput',false)), '0'];

pname = {'surfaceFluxRampFactor'};
plistarr = ['1'];

% punit = '';
useparam = 'on';
pcontinuationmode = 'no';

          
studyNum = m.m.study.size() + 1;
study_id = sprintf('study%u',studyNum);

m.m.study.create(study_id);

solNum = m.m.sol.size() + 1;
sol_id = sprintf('sol%u',solNum);

m.m.sol.create(sol_id);
m.m.sol(sol_id).study(study_id); % ?
m.m.sol(sol_id).attach(study_id); % ?

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

% set scaling for variables
m.m.sol(sol_id).feature(variables_id).set('scalemethod','init');

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

m.m.sol(sol_id).feature(solver_id).set('probesel', 'none');
m.m.sol(sol_id).feature(solver_id).set('nonlin', 'on');
m.m.sol(sol_id).feature(solver_id).set('stol', '1e-3');

m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'const');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermconst', 'itertol');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '30');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname); 

display('Setting individual active physics for step...');

m.m.study(study_id).feature(studyStep_id).set('disabledphysics', {''});
m.m.study(study_id).feature(studyStep_id).activate('NernstPlanckEquation', false);
m.m.study(study_id).feature(studyStep_id).activate('PoissonEquation', false);
% m.m.study(study_id).feature(studyStep_id).activate('WeakFormulation1d', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroSurfaceCurrent', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroNetCurrent0d', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroNetCurrent1d', false);

m.saveState;
%% solve 1d 
m.m.sol(sol_id).run(store_id);

%% export results of 1d solution and import as interpolated function
storedSolutions = m.m.study(study_id).getSolverSequences('Stored');
solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(1)));
dset = char(solInfo.dataset(2));

m.updateFluxEvaluations1d(dset);

txtFileName = m.exportSolution1d(dset);
m.m.func('interpSolution1d').set('filename', txtFileName);
m.m.func('interpSolution1d').set('nargs', '1');
m.m.func('interpSolution1d').active(true);

return;

%% plot 1d results
% m.updateEvaluations1d(dset);
m.plot1dSolution(dset);

% return;

%% meshing for 2d component

% mapped mesh
m.m.disableUpdates(false); % necessary before meshing
m.hMaxFactor = 0.1;
m.mesh1D();
% m.createMappedMeshForSimpleGeometry();
% m.updateMappedMeshForSimpleGeometry();
m.createMappedMesh();
m.updateMappedMesh();
m.m.mesh('mappedMesh').run;
m.saveState;

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

% set scaling for variables
m.m.sol(sol_id).feature(variables_id).set('scalemethod','init');

% what do those stand for:
m.m.study(study_id).feature(studyStep_id).set('initstudyhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('initsolhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('solnumhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('notstudyhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('notsolhide', 'on');
m.m.study(study_id).feature(studyStep_id).set('notsolnumhide', 'on');

m.m.study(study_id).feature(studyStep_id).set('mesh', {'geom' 'mappedMesh' 'geom1d' 'nomesh'});

m.m.study(study_id).feature(studyStep_id).set('useparam', useparam);
m.m.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
m.m.study(study_id).feature(studyStep_id).set('pname', pname);
m.m.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
m.m.study(study_id).feature(studyStep_id).set('pcontinuationmode', pcontinuationmode);
m.m.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
m.m.study(study_id).feature(studyStep_id).set('showdistribute', true);


m.m.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);
m.m.sol(sol_id).feature(solver_id).set('probesel', 'none');
m.m.sol(sol_id).feature(solver_id).set('nonlin', 'on');
m.m.sol(sol_id).feature(solver_id).set('stol', '1e-3');

% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'const');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermconst', 'itertol');
m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '15');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
% m.m.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
m.m.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname); 

% display('Setting individual active physics for step...');

m.m.study(study_id).feature(studyStep_id).set('disabledphysics', {''});
% m.m.study(study_id).feature(studyStep_id).activate('WeakFormulation', false);
m.m.study(study_id).feature(studyStep_id).activate('WeakFormulation1d', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroSurfaceCurrent', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroNetCurrent0d', false);
m.m.study(study_id).feature(studyStep_id).activate('zeroNetCurrent1d', false);

m.saveState;

return

%% plot single-step 2d
storedSolutions = m.m.study(study_id).getSolverSequences('Stored');
solInfo = mphsolinfo(m.m,'soltag',char(storedSolutions(2)));
dset = char(solInfo.dataset(1));

m.plotPresetCrossections(dset);
m.sweepHorizontalCrossection(dset,0.05)
m.sweepVerticalCrossection(dset,m.w/10,0.3);
m.update2dPlots(dset);

m.update2dParametricPlots(dset,'deltaPhi');

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

return;

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

tag = '2016_06_15_21_32_20_duval2001bipolar';

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
