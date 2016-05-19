m = jlh.BpeModel;
m.newProject('duval2001bipolar_case');
% m.setSampleCaseParameters();
m.setDuval2001BipolarCaseParameters(1);
m.prepareIdentifiers;
% E_m = m.calcMixedPotential;
m.PHI_bpe = 0;
% fprintf('Mixed potential: %f',E_m);

%% m.createModel;
m.m = ModelUtil.create(m.model_tag);   % creates model on COMSOL server

m.createFunctions()
m.m.modelNode.create(m.comp_id); 
m.updateParameters();
m.makeChoppedGeometry();

%% 
% m.createGeometry();
% m.createSelections();
m.createOperators()
m.createVariables();
m.createPhysics();
m.createMesh();
m.createProbes();
m.create1dPlots();
m.create2dPlots();
m.standardView = m.m.view.create('standardView', 'geom');
m.ddlView      = m.m.view.create('ddlView', 'geom');
m.createDatasets();

%% m.updateModel;
m.updateFunctions();
% m.updateSelections();
% m.updateOperators();
m.updateOperatorsForChoppedGeometry();
m.updatePhysicsForChoppedGeometry();
%%

m.hMaxFactor = 1;
m.updateChoppedMesh();

% m.updateMeshMapped();
% m.addParametricSweepStudy('chargeDensityRampFactor',0);
% m.addParametricSweepStudy({'chargeDensityRampFactor','phiSymmetryFactor'},{[0,1]; [0:0.1:1]});
m.addParametricSweepStudy({'deltaPhiRampFactor'},{1});
% m.setUniformFeederElectrodePotentialBC
% m.setUniformFeederElectrodePotentialAndLinearBulkPotentialBC
% m.setFeederElectrodeCurrentFluxAndBPEPotentialBC
m.setMinimalBC;
m.setUniformFeederElectrodePotentialAndLinearBulkPotentialBC;
m.activateBpe;
%% meshing
m.testingMesh.run
m.saveState

%% run
m.stationaryStudy.run

%% after getting initial values or result
dset = m.getLatestDataset();

m.updateDatasets(dset)
% m.iterateStandardPlots
m.update1dPlots
m.updateSurfacePlots(dset)
%% save 2d plots
m.updateViews
m.update2dPlots(dset)

%% parametric plots
m.updateGlobalPlotsBySolnum(dset,'phiSymmetryFactor');
m.update1dParametricPlots(dset);
m.update2dParametricPlots(dset);