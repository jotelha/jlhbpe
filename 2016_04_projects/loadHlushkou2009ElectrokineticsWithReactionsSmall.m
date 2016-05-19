m = jlh.BpeModel;
% m.newProject('Dhopeshwarkar2008ElectrokineticsCaseSmall');
% m.setDhopeshwarkar2008ElectrokineticsCaseParameters;

setHlushkou2009ElectrokineticsWithReactionsParametersSmall

m.prepareIdentifiers;

%% load;

% listing = dir('dat');
% tag = listing(end).name;

tag = '2016_04_21_13_40_39_Hlushkou2009ElectrokineticsWithReactionsParametersSmall';

m.projectName = tag;
m.projectPath = ['dat\',tag];

tags = ModelUtil.tags;
m.loadFromServer( char( tags(1) ));

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
m.updateGlobalPlotsBySolnum(dset,'phi_bpe_init');
m.update1dParametricPlots(dset);
m.update2dParametricPlots(dset);