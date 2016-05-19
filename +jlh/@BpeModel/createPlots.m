function createPlots(obj)
     %% create standard plots, moved here on 2016/01/27
    obj.plotExporter1d      = obj.m.result.export.create('plotExporter1d', 'Image1D');
    
    %% create potential plots
    fprintf('  Creating empty plots...\n');

    obj.standardPhiPlotGroup = obj.m.result.create('standardPhiPlotGroup',1);
    obj.standardPhiPlot = obj.standardPhiPlotGroup.create('standardPhiPlot','LineGraph');
    obj.phiPBPlot = obj.standardPhiPlotGroup.feature.duplicate('phiPBPlot','standardPhiPlot');


    %% create concentrations plot
    obj.standardConcentrationsPlotGroup = obj.m.result.create('standardConcentrationsPlotGroup',1);
    obj.standardConcentrationsPlot = cell(1,obj.numberOfSpecies);
    for i=1:obj.numberOfSpecies
        obj.standardConcentrationsPlot{i} = obj.standardConcentrationsPlotGroup.create(obj.standardConcentrationsPlot_id{i},'LineGraph');
    end

    %% create log plots
    obj.logConcentrationsPlotGroup = obj.m.result.create('logConcentrationsPlotGroup',1);
    obj.logConcentrationsPlot = cell(1,obj.numberOfSpecies);
    for i=1:obj.numberOfSpecies
        obj.logConcentrationsPlot{i} = obj.logConcentrationsPlotGroup.create(obj.logConcentrationsPlot_id{i},'LineGraph');
    end

%     %% create plots of ddl region just at surface, phi
%     % potential
%     obj.phiAtSurfacePlotGroup = obj.m.result.create('phiAtSurfacePlotGroup',1);
%     obj.phiAtSurfacePlot = obj.phiAtSurfacePlotGroup.create('phiAtSurfacePlot','LineGraph');
%     obj.phiAtSurfacePBPlot = obj.phiAtSurfacePlotGroup.feature.duplicate('phiAtSurfacePBPlot','phiAtSurfacePlot');
% 
%     %% create plots of ddl region just at surface, log concentrations
%     obj.logConcentrationsAtSurfacePlotGroup = obj.m.result.create('logConcentrationsAtSurfacePlotGroup',1);
%     obj.logConcentrationsAtSurfacePlot = cell(1,obj.numberOfSpecies);
%     obj.logConcentrationsAtSurfacePlot_id = nprintf('logConcentrationsAtSurfacePlot%d',obj.numberOfSpecies);
%     obj.logPBConcentrationsAtSurfacePlot = cell(1,obj.numberOfSpecies);
%     obj.logPBConcentrationsAtSurfacePlot_id = nprintf('logPBConcentrationsAtSurfacePlot%d',obj.numberOfSpecies);
%     for i=1:obj.numberOfSpecies
%         obj.logConcentrationsAtSurfacePlot{i} = obj.logConcentrationsAtSurfacePlotGroup.create(obj.logConcentrationsAtSurfacePlot_id{i},'LineGraph');
%         obj.logPBConcentrationsAtSurfacePlot{i} = obj.logConcentrationsAtSurfacePlotGroup.create(obj.logPBConcentrationsAtSurfacePlot_id{i},'LineGraph');
%     end
% 
%     fprintf('  Finished creating COMSOL model.\n');

    %% create 2d plots
    obj.create2dPlots();
    
    %% create views
   obj.standardView = obj.m.view.create('standardView', 'geom');
   obj.ddlView      = obj.m.view.create('ddlView', 'geom');
end