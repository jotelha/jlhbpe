m = jlh.BpeModel;
% m.newProject('sample_project');
m.setSampleCaseParameters();
m.prepareIdentifiers;

m.projectName = '2016_02_28_15_15_34_massive_sweep';
m.loadMph;
% E_m = m.calcMixedPotential;
% m.PHI_bpe = E_m;
% % fprintf('Mixed potential: %f',E_m);
% m.createModel;
% m.updateModel;
% % m.addParametricSweepStudy('chargeDensityRampFactor',0);
% m.addParametricSweepStudy({'chargeDensityRampFactor','phiSymmetryFactor'},{[0,1]; [0:0.1:1]});
% % m.setUniformFeederElectrodePotentialBC
% % m.setUniformFeederElectrodePotentialAndLinearBulkPotentialBC
% % m.setFeederElectrodeCurrentFluxAndBPEPotentialBC
% m.setMinimalBC;
% m.activateBpe;
% %% meshing
% m.testingMesh.run
% m.saveState
% 
% %% run
% m.stationaryStudy.run

%% afte getting initial values or result
dset = m.getLatestDataset();

m.updateDatasets(dset)
% m.iterateStandardPlots
m.update1dPlots

%% save 2d plots
m.updateViews
m.update2dPlots(dset)

%% parametric plots
m.updateGlobalPlotsBySolnum(dset,'phiSymmetryFactor');
m.update1dParametricPlots(dset);
m.update2dParametricPlots(dset);