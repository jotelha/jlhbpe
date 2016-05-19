m = jlh.BpeModel;
m.newProject('Hlushkou2009ElectrokineticsWithReactionsParametersSmall');

% server log file
logFile = [pwd(),'\',m.projectPath,'\comsol.log'];
% fclose(fopen(logFile, 'w'));
ModelUtil.showProgress(logFile);
% m.setDhopeshwarkar2008ElectrokineticsCaseParameters;

setHlushkou2009ElectrokineticsWithReactionsParametersSmall2

E_m = m.calcMixedPotential;
m.plotCurrents;
m.PHI_bpe = E_m;

m.prepareIdentifiers;

%% m.createModel;
loadedModels    = ModelUtil.tags;
isLoaded        = arrayfun( @(s) strcmp(m.model_tag,s),loadedModels);
if any(isLoaded)
    ModelUtil.remove(m.model_tag);
end
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
% 
m.hMaxFactor = 0.01;
m.updateChoppedMesh();
m.replicateMeshPrototype;
m.finalizeChoppedMesh;
m.standardMesh.run;
m.saveState;

%% run
m.stationaryStudy.run

%% a second study?
m.addBpeStudy({'deltaPhiRampFactor'},{1});
phi_bpe = -30:1:70;
m.addBpeStudy({'phi_bpe'},{phi_bpe});


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