m = jlh.BpeModel;
m.newProject('duval2001bipolar_case');
m.setDuval2001BipolarCaseParameters(1);
m.prepareIdentifiers;
E_m = m.calcMixedPotential;
%m.PHI_bpe = E_m;
%%
% fprintf('Mixed potential: %f',E_m);
m.create1dModel;
m.hMaxFactor = 1e-3;
m.update1dModel;
% m.addParametricSweepStudy('chargeDensityRampFactor',0);
par = cell(1,1);
par{1} = ((0.9*m.phi_bpe):(m.phi_bpe/100):(1.1*m.phi_bpe));
m.addParametricSweepStudy({'phi_bpe'},par);
% m.setUniformFeederElectrodePotentialBC
% m.setUniformFeederElectrodePotentialAndLinearBulkPotentialBC
% m.setFeederElectrodeCurrentFluxAndBPEPotentialBC
m.setMinimalBC;
m.activateBpe;
%% meshing
m.standardMesh.run
m.saveState

%% run
m.stationaryStudy.run

%% afte getting initial values or result
dset = m.getLatestDataset();
m.updateGlobalPlotsBySolnum(dset,'phi_bpe');

% m.updateDatasets(dset)
% % m.iterateStandardPlots
% m.update1dPlots
% 
% %% save 2d plots
% m.updateViews
% m.update2dPlots(dset)
% 
% %% parametric plots
% m.updateGlobalPlots;
% m.update1dParametricPlots(dset);
% m.update2dParametricPlots(dset);