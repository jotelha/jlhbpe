m = jlh.BpeModel;
m.newProject('duval2001bipolar_case');
% m.setSampleCaseParameters();
m.setDuval2001BipolarCaseParameters(1);
m.prepareIdentifiers;
% E_m = m.calcMixedPotential;
m.PHI_bpe = 0;
% fprintf('Mixed potential: %f',E_m);
m.createModel;
%% m.updateModel;
m.updateParameters();
m.updateFunctions();
m.updateGeometry();
m.updateSelections();
m.updateOperators();
m.updatePhysics();

m.hMaxFactor = 0.1;
m.updateMesh();

% m.updateMeshMapped();
% m.addParametricSweepStudy('chargeDensityRampFactor',0);
% m.addParametricSweepStudy({'chargeDensityRampFactor','phiSymmetryFactor'},{[0,1]; [0:0.1:1]});
m.addParametricSweepStudy({'deltaPhiRampFactor'},{[1:1:3]});
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

%% save 2d plots
m.updateViews
m.update2dPlots(dset)

%% parametric plots
m.updateGlobalPlotsBySolnum(dset,'phiSymmetryFactor');
m.update1dParametricPlots(dset);
m.update2dParametricPlots(dset);