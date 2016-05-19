function obj = create2dPlots(obj)
    obj.plotExporter2d      = obj.m.result.export.create('plotExporter2d', 'Image2D');
    
    obj.phiSurfacePlotGroup = obj.m.result.create('phiSurfacePlotGroup', 'PlotGroup2D');
    obj.cSurfacePlotGroup   = obj.m.result.create('cSurfacePlotGroup', 'PlotGroup2D');
    
    obj.phiContourPlotGroup = obj.m.result.create('phiContourPlotGroup', 'PlotGroup2D');
    obj.cContourPlotGroup   = obj.m.result.create('cContourPlotGroup', 'PlotGroup2D');
   
    obj.phiStreamlinePlotGroup  = obj.m.result.create('phiStreamlinePlotGroup', 'PlotGroup2D');
    obj.cStreamlinePlotGroup    = obj.m.result.create('cStreamlinePlotGroup', 'PlotGroup2D');

    obj.phiArrowSurfacePlotGroup= obj.m.result.create('phiArrowSurfacePlotGroup', 'PlotGroup2D');
    obj.cArrowSurfacePlotGroup  = obj.m.result.create('cArrowSurfacePlotGroup', 'PlotGroup2D');
    
    obj.phiArrowLinePlotGroup   = obj.m.result.create('phiArrowLinePlotGroup', 'PlotGroup2D');
    obj.cArrowLinePlotGroup     = obj.m.result.create('cArrowLinePlotGroup', 'PlotGroup2D');
    
    obj.meshPlotGroup           = obj.m.result.create('meshPlotGroup', 'PlotGroup2D');
%     obj.m.result.create('pg9', 'PlotGroup2D');
  
    obj.m.result('phiSurfacePlotGroup').create('phiSurfacePlot', 'Surface');
    obj.m.result('cSurfacePlotGroup').create('cSurfacePlot', 'Surface');
    
    obj.m.result('phiContourPlotGroup').create('phiContourPlot', 'Contour');
    obj.m.result('cContourPlotGroup').create('cContourPlot', 'Contour');
    
    obj.m.result('phiStreamlinePlotGroup').create('phiStreamlinePlot', 'Streamline');
    obj.m.result('cStreamlinePlotGroup').create('cStreamlinePlot', 'Streamline');
    % obj.m.result('phiStreamlinePlotGroup').feature('phiStreamlinePlot').selection.all;
    
    obj.m.result('phiArrowSurfacePlotGroup').create('phiArrowSurfacePlot', 'ArrowSurface');
    obj.m.result('cArrowSurfacePlotGroup').create('cArrowSurfacePlot', 'ArrowSurface');
    
    obj.m.result('phiArrowLinePlotGroup').create('phiArrowLinePlot', 'ArrowLine');
    obj.m.result('cArrowLinePlotGroup').create('cArrowLinePlot', 'ArrowLine');

    obj.m.result('meshPlotGroup').create('meshPlot', 'Mesh');
    
    % conductivity
    obj.kappaSurfacePlotGroup = obj.m.result.create('kappaSurfacePlotGroup', 'PlotGroup2D');
    obj.kappaSurfacePlotGroup.create('kappaSurfacePlot', 'Surface');
    obj.kappaContourPlotGroup = obj.m.result.create('kappaContourPlotGroup', 'PlotGroup2D');
    obj.kappaContourPlotGroup.create('kappaContourPlot', 'Contour');

    % species flux plots
    obj.nStreamlinePlotGroup    = obj.m.result.create('nStreamlinePlotGroup', 'PlotGroup2D');
    obj.nArrowSurfacePlotGroup  = obj.m.result.create('nArrowSurfacePlotGroup', 'PlotGroup2D');
    obj.nArrowLinePlotGroup     = obj.m.result.create('nArrowLinePlotGroup', 'PlotGroup2D'); % species flux plots
    
    obj.nStreamlinePlotGroup.create('nStreamlinePlot', 'Streamline');
    obj.nArrowSurfacePlotGroup.create('nArrowSurfacePlot', 'ArrowSurface');
    obj.nArrowLinePlotGroup.create('nArrowLinePlot', 'ArrowLine');
    
    % current density plots
    obj.iStreamlinePlotGroup    = obj.m.result.create('iStreamlinePlotGroup', 'PlotGroup2D');
    obj.iArrowSurfacePlotGroup  = obj.m.result.create('iArrowSurfacePlotGroup', 'PlotGroup2D');
    obj.iArrowLinePlotGroup     = obj.m.result.create('iArrowLinePlotGroup', 'PlotGroup2D');
    
    obj.iStreamlinePlotGroup.create('iStreamlinePlot', 'Streamline');
    obj.iArrowSurfacePlotGroup.create('iArrowSurfacePlot', 'ArrowSurface');
    obj.iArrowLinePlotGroup.create('iArrowLinePlot', 'ArrowLine');
%     obj.m.result.export.create('img1', 'Image2D');
%     obj.m.result.export.create('anim1', 'Animation');
end