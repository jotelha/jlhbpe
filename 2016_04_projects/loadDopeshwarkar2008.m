m = jlh.BpeModel;
% m.newProject('Dhopeshwarkar2008ElectrokineticsCaseSmall');
% m.setDhopeshwarkar2008ElectrokineticsCaseParameters;

setDopeshwarkar2008Electrokinetics

m.prepareIdentifiers;

%% load;

% listing = dir('dat');
% tag = listing(end).name;
tag = '2016_04_13_16_04_54_Dhopeshwarkar2008ElectrokineticsCase';

m.projectName = tag;
m.projectPath = ['dat\',tag];

m.loadFromServer(tag);

%% complte the mesh
m.replicateTestingMeshPrototype;
m.finalizeChoppedTestingMesh;

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