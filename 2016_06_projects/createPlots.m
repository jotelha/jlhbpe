%% plot exports
model.result.export.create('plotExporter1d', 'Image1D');
model.result.export.create('plotExporter2d', 'Image2D');

%% 1d plots
% create potential plots
fprintf('  Creating empty plots...\n');

model.result.create('standardPhiPlotGroup',1);
model.result('standardPhiPlotGroup').create('standardPhiPlot','LineGraph');
%model.result('standardPhiPlotGroup').feature.duplicate('phiPBPlot','standardPhiPlot');


% create concentrations plot
model.result.create('standardConcentrationsPlotGroup',1);
for i=1:obj.numberOfSpecies
	model.result('standardConcentrationsPlotGroup').create(obj.standardConcentrationsPlot_id{i},'LineGraph');
end

% create log plots
model.result.create('logConcentrationsPlotGroup',1);
for i=1:obj.numberOfSpecies
	model.result('logConcentrationsPlotGroup').create(obj.logConcentrationsPlot_id{i},'LineGraph');
end

model.result.create('globalPlotGroup','PlotGroup1D');
obj.globalPlotGroup.create('globalPlot','Global');
	
model.result.create('multiPurpose1dPlotGroup','PlotGroup1D');
for j=1:obj.nMultiPurpose1dPlots
	model.result('multiPurpose1dPlotGroup').create(obj.multiPurpose1dPlot_id{j}, 'LineGraph');
end
	
%% 2d plots
model.result.create('phiSurfacePlotGroup', 'PlotGroup2D');
model.result.create('cSurfacePlotGroup', 'PlotGroup2D');

model.result.create('phiContourPlotGroup', 'PlotGroup2D');
model.result.create('cContourPlotGroup', 'PlotGroup2D');

model.result.create('phiStreamlinePlotGroup', 'PlotGroup2D');
model.result.create('cStreamlinePlotGroup', 'PlotGroup2D');

model.result.create('phiArrowSurfacePlotGroup', 'PlotGroup2D');
model.result.create('cArrowSurfacePlotGroup', 'PlotGroup2D');

model.result.create('phiArrowLinePlotGroup', 'PlotGroup2D');
model.result.create('cArrowLinePlotGroup', 'PlotGroup2D');

model.result.create('meshPlotGroup', 'PlotGroup2D');

model.result('phiSurfacePlotGroup').create('phiSurfacePlot', 'Surface');
model.result('cSurfacePlotGroup').create('cSurfacePlot', 'Surface');

model.result('phiContourPlotGroup').create('phiContourPlot', 'Contour');
model.result('cContourPlotGroup').create('cContourPlot', 'Contour');

model.result('phiStreamlinePlotGroup').create('phiStreamlinePlot', 'Streamline');
model.result('cStreamlinePlotGroup').create('cStreamlinePlot', 'Streamline');

model.result('phiArrowSurfacePlotGroup').create('phiArrowSurfacePlot', 'ArrowSurface');
model.result('cArrowSurfacePlotGroup').create('cArrowSurfacePlot', 'ArrowSurface');

model.result('phiArrowLinePlotGroup').create('phiArrowLinePlot', 'ArrowLine');
model.result('cArrowLinePlotGroup').create('cArrowLinePlot', 'ArrowLine');

model.result('meshPlotGroup').create('meshPlot', 'Mesh');

% conductivity
model.result.create('kappaSurfacePlotGroup', 'PlotGroup2D');
model.result('kappaSurfacePlotGroup').create('kappaSurfacePlot', 'Surface');
model.result('kappaSurfacePlotGroup').create('kappaContourPlotGroup', 'PlotGroup2D');
model.result('kappaSurfacePlotGroup').create('kappaContourPlot', 'Contour');

% species flux plots
model.result.create('nStreamlinePlotGroup', 'PlotGroup2D');
model.result.create('nArrowSurfacePlotGroup', 'PlotGroup2D');
model.result.create('nArrowLinePlotGroup', 'PlotGroup2D'); % species flux plots

model.result('nStreamlinePlotGroup').create('nStreamlinePlot', 'Streamline');
model.result('nArrowSurfacePlotGroup').create('nArrowSurfacePlot', 'ArrowSurface');
model.result('nArrowLinePlotGroup').create('nArrowLinePlot', 'ArrowLine');

% current density plots
model.result.create('iStreamlinePlotGroup', 'PlotGroup2D');
model.result.create('iArrowSurfacePlotGroup', 'PlotGroup2D');
model.result.create('iArrowLinePlotGroup', 'PlotGroup2D');

model.result('iStreamlinePlotGroup').create('iStreamlinePlot', 'Streamline');
model.result('iArrowSurfacePlotGroup.')create('iArrowSurfacePlot', 'ArrowSurface');
model.result('iArrowLinePlotGroup').create('iArrowLinePlot', 'ArrowLine');