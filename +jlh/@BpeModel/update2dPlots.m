function obj = update2dPlots(obj,dset)
    % try to use current path as root, might fail on server?
    fullProjectPath = [pwd(),'\',obj.projectPath];
    
    solSubFolder = [fullProjectPath,'\',dset];
    if( ~exist( solSubFolder,'dir') )
        mkdir(solSubFolder);
    end
    
    standard2dSubFolder = [solSubFolder,'\standard2d'];
    if( ~exist( standard2dSubFolder,'dir') )
        mkdir(standard2dSubFolder);
    end
    fileNameTemplate = [strrep(standard2dSubFolder,'\','\\'),'\\%s.png'];

    streamlinePosmeth = 'uniform'; % magnitude
    streamlineDist = 0.02;
    arrowScale = 0.01;
    arrowScaleActive = false;
%     arrowScaleActive = true;
%     arrowLength = 'normalized'; % 'logarithmic'
    arrowLength = 'proportional'; % 'logarithmic'

%% setup exporter
%     obj.m.result.export('plotExporter2d').set('plotgroup', 'pg3');
%     obj.m.result.export('plotExporter2d').set('pngfilename', 'G:\johnny\matlab\2016_02_projects\dat\2016_02_17_17_06_56_sample_project\exportTest.png');
    obj.m.result.export('plotExporter2d').set('view', 'standardView');
    obj.m.result.export('plotExporter2d').set('printunit', 'mm');
    obj.m.result.export('plotExporter2d').set('webunit', 'px');
    obj.m.result.export('plotExporter2d').set('printheight', '90');
    obj.m.result.export('plotExporter2d').set('webheight', '600');
    obj.m.result.export('plotExporter2d').set('printwidth', '120');
    obj.m.result.export('plotExporter2d').set('webwidth', '800');
    obj.m.result.export('plotExporter2d').set('printlockratio', 'off');
    obj.m.result.export('plotExporter2d').set('weblockratio', 'off');
    obj.m.result.export('plotExporter2d').set('printresolution', '300');
    obj.m.result.export('plotExporter2d').set('webresolution', '96');
    obj.m.result.export('plotExporter2d').set('size', 'manualprint');
    obj.m.result.export('plotExporter2d').set('antialias', 'on');
    obj.m.result.export('plotExporter2d').set('zoomextents', 'off'); % off
    obj.m.result.export('plotExporter2d').set('title', 'on');
    obj.m.result.export('plotExporter2d').set('legend', 'on');
    obj.m.result.export('plotExporter2d').set('logo', 'on');
    obj.m.result.export('plotExporter2d').set('options', 'on');
    obj.m.result.export('plotExporter2d').set('fontsize', '9');
    obj.m.result.export('plotExporter2d').set('customcolor', [1 1 1]);
    obj.m.result.export('plotExporter2d').set('background', 'color');
    obj.m.result.export('plotExporter2d').set('axes', 'on');
    obj.m.result.export('plotExporter2d').set('qualitylevel', '92');
    obj.m.result.export('plotExporter2d').set('qualityactive', 'on');
    obj.m.result.export('plotExporter2d').set('imagetype', 'png');
    
%% potential plots
    obj.m.result('phiSurfacePlotGroup').set('data',dset);
    obj.m.result('phiSurfacePlotGroup').set('view', 'standardView');
    obj.m.result('phiSurfacePlotGroup').feature('phiSurfacePlot').label('phi');
    obj.m.result('phiSurfacePlotGroup').feature('phiSurfacePlot').set('descr', 'Dependent variable phi');
    obj.m.result('phiSurfacePlotGroup').feature('phiSurfacePlot').set('expr', 'phi');
    obj.m.result('phiSurfacePlotGroup').run;

    obj.m.result.export('plotExporter2d').set('plotgroup', 'phiSurfacePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'phiSurfacePlotGroup'));
    obj.m.result.export('plotExporter2d').run();
    
    obj.m.result('phiContourPlotGroup').set('data',dset);
    obj.m.result('phiContourPlotGroup').set('view', 'standardView');
    obj.m.result('phiContourPlotGroup').feature('phiContourPlot').label('phi');
    obj.m.result('phiContourPlotGroup').feature('phiContourPlot').set('descr', 'Dependent variable phi');
    obj.m.result('phiContourPlotGroup').feature('phiContourPlot').set('expr', 'phi');
    obj.m.result.export('plotExporter2d').set('plotgroup', 'phiContourPlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'phiContourPlotGroup'));
    obj.m.result.export('plotExporter2d').run();

    obj.m.result('phiStreamlinePlotGroup').set('data',dset);
    obj.m.result('phiStreamlinePlotGroup').set('view', 'standardView');
    obj.m.result('phiStreamlinePlotGroup').feature('phiStreamlinePlot').set('descr', 'Gradient of phi');
    obj.m.result('phiStreamlinePlotGroup').feature('phiStreamlinePlot').set('posmethod', streamlinePosmeth);
    obj.m.result('phiStreamlinePlotGroup').feature('phiStreamlinePlot').set('udist', streamlineDist);
%     obj.m.result('phiStreamlinePlotGroup').feature('phiStreamlinePlot').set('mdensity', '15');
    obj.m.result('phiStreamlinePlotGroup').feature('phiStreamlinePlot').set('expr', {'phix' 'phiy'});
    obj.m.result.export('plotExporter2d').set('plotgroup', 'phiStreamlinePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'phiStreamlinePlotGroup'));
    obj.m.result.export('plotExporter2d').run();
   
    obj.m.result('phiArrowSurfacePlotGroup').set('data',dset);
    obj.m.result('phiArrowSurfacePlotGroup').set('view', 'standardView');
    obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('xnumber', '20');
   	obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('ynumber', '5');
   	obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('expr', {'phix' 'phiy'});
    obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('arrowlength', arrowLength);
   	obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('scale', arrowScale);
   	obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('descr', 'Gradient of phi');
   	obj.m.result('phiArrowSurfacePlotGroup').feature('phiArrowSurfacePlot').set('scaleactive', arrowScaleActive);
    obj.m.result.export('plotExporter2d').set('plotgroup', 'phiArrowSurfacePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'phiArrowSurfacePlotGroup'));
    obj.m.result.export('plotExporter2d').run();
    
    obj.m.result('phiArrowLinePlotGroup').set('data',dset);
    obj.m.result('phiArrowLinePlotGroup').set('view', 'standardView');
    obj.m.result('phiArrowLinePlotGroup').feature('phiArrowLinePlot').set('expr', {'phix' 'phiy'});
    obj.m.result('phiArrowLinePlotGroup').feature('phiArrowLinePlot').set('arrowlength', arrowLength);
    obj.m.result('phiArrowLinePlotGroup').feature('phiArrowLinePlot').set('scale', arrowScale);
    obj.m.result('phiArrowLinePlotGroup').feature('phiArrowLinePlot').set('descr', 'Gradient of phi');
    obj.m.result('phiArrowLinePlotGroup').feature('phiArrowLinePlot').set('scaleactive', arrowScaleActive);
%     obj.m.result('obj.phiArrowLinePlotGroup').feature('phiArrowLinePlot').set('wireframecolor', 'none');
    obj.m.result.export('plotExporter2d').set('plotgroup', 'phiArrowLinePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'phiArrowLinePlotGroup'));
    obj.m.result.export('plotExporter2d').run();
   
    %% concentration plots
    for i=1:obj.numberOfSpecies
        
        % species flux N    
        fileNameTemplate = [strrep(standard2dSubFolder,'\','\\'),'\\%s_',obj.c_id{i},'.png'];
%         fileNameTemplate = [fullProjectPath,'\\%s_',obj.c_id{i},'.png'];

        obj.m.result('nStreamlinePlotGroup').set('data',dset);
        obj.m.result('nStreamlinePlotGroup').set('view', 'standardView');
        obj.m.result('nStreamlinePlotGroup').feature('nStreamlinePlot').set('descr', 'species flux');
        obj.m.result('nStreamlinePlotGroup').feature('nStreamlinePlot').set('posmethod', streamlinePosmeth);
        obj.m.result('nStreamlinePlotGroup').feature('nStreamlinePlot').set('udist', streamlineDist);
%         obj.m.result('nStreamlinePlotGroup').feature('nStreamlinePlot').set('mdensity', '10');
        obj.m.result('nStreamlinePlotGroup').feature('nStreamlinePlot').set('expr', {obj.nx_id{i} obj.ny_id{i}});
        obj.m.result.export('plotExporter2d').set('plotgroup', 'nStreamlinePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'nStreamlinePlotGroup'));
        obj.m.result.export('plotExporter2d').run();
        
        obj.m.result('nArrowSurfacePlotGroup').set('data',dset);
        obj.m.result('nArrowSurfacePlotGroup').set('view', 'standardView');
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('xnumber', '20');
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('ynumber', '5');
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('expr', {obj.nx_id{i} obj.ny_id{i}});
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('arrowlength', arrowLength);
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('scale', arrowScale);
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('descr', 'species flux');
        obj.m.result('nArrowSurfacePlotGroup').feature('nArrowSurfacePlot').set('scaleactive', arrowScaleActive);
        obj.m.result.export('plotExporter2d').set('plotgroup', 'nArrowSurfacePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'nArrowSurfacePlotGroup'));
        obj.m.result.export('plotExporter2d').run();
        
        obj.m.result('nArrowLinePlotGroup').set('data',dset);
        obj.m.result('nArrowLinePlotGroup').set('view', 'standardView');
        obj.m.result('nArrowLinePlotGroup').feature('nArrowLinePlot').set('expr', {obj.nx_id{i} obj.ny_id{i}});
        obj.m.result('nArrowLinePlotGroup').feature('nArrowLinePlot').set('arrowlength', arrowLength);
        obj.m.result('nArrowLinePlotGroup').feature('nArrowLinePlot').set('scale', arrowScale);
        obj.m.result('nArrowLinePlotGroup').feature('nArrowLinePlot').set('descr', 'Gradient of c');
        obj.m.result('nArrowLinePlotGroup').feature('nArrowLinePlot').set('scaleactive', arrowScaleActive);
        obj.m.result.export('plotExporter2d').set('plotgroup', 'nArrowLinePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'nArrowLinePlotGroup'));
        obj.m.result.export('plotExporter2d').run();

        
        % species concentration c
%         fileNameTemplate = [fullProjectPath,'\\%s_',obj.c_id{i},'.png'];

        
        obj.m.result('cSurfacePlotGroup').set('data',dset);
        obj.m.result('cSurfacePlotGroup').set('view', 'standardView');
        obj.m.result('cSurfacePlotGroup').feature('cSurfacePlot').label(obj.c_id{i});
        obj.m.result('cSurfacePlotGroup').feature('cSurfacePlot').set('descr', obj.c_id{i});
        obj.m.result('cSurfacePlotGroup').feature('cSurfacePlot').set('expr', obj.c_id{i});
        obj.m.result.export('plotExporter2d').set('plotgroup', 'cSurfacePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'cSurfacePlotGroup'));
        obj.m.result.export('plotExporter2d').run();

    %     model.result('pg7').feature('con1').set('unit', '');
        obj.m.result('cContourPlotGroup').set('data',dset);
        obj.m.result('cContourPlotGroup').set('view', 'standardView');
        obj.m.result('cContourPlotGroup').feature('cContourPlot').label(sprintf('log(%s)',obj.c_id{i}));
        obj.m.result('cContourPlotGroup').feature('cContourPlot').set('descr', sprintf('log(%s)',obj.c_id{i}));
        obj.m.result('cContourPlotGroup').feature('cContourPlot').set('expr', sprintf('log(%s)',obj.c_id{i}));
        obj.m.result.export('plotExporter2d').set('plotgroup', 'cContourPlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'cContourPlotGroup'));
        obj.m.result.export('plotExporter2d').run();

        obj.m.result('cStreamlinePlotGroup').set('data',dset);
        obj.m.result('cStreamlinePlotGroup').set('view', 'standardView');
        obj.m.result('cStreamlinePlotGroup').feature('cStreamlinePlot').set('descr', 'Gradient of c');
        obj.m.result('cStreamlinePlotGroup').feature('cStreamlinePlot').set('posmethod', streamlinePosmeth);
        obj.m.result('cStreamlinePlotGroup').feature('cStreamlinePlot').set('udist', streamlineDist);
%         obj.m.result('cStreamlinePlotGroup').feature('cStreamlinePlot').set('mdensity', '10');
        obj.m.result('cStreamlinePlotGroup').feature('cStreamlinePlot').set('expr', {obj.cx_id{i} obj.cy_id{i}});
        obj.m.result.export('plotExporter2d').set('plotgroup', 'cStreamlinePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'cStreamlinePlotGroup'));
        obj.m.result.export('plotExporter2d').run();

        obj.m.result('cArrowSurfacePlotGroup').set('data',dset);
        obj.m.result('cArrowSurfacePlotGroup').set('view', 'standardView');
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('xnumber', '20');
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('ynumber', '5');
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('expr', {obj.cx_id{i} obj.cy_id{i}});
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('arrowlength', arrowLength);
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('scale', arrowScale);
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('descr', 'Gradient of c');
        obj.m.result('cArrowSurfacePlotGroup').feature('cArrowSurfacePlot').set('scaleactive', arrowScaleActive);
        obj.m.result.export('plotExporter2d').set('plotgroup', 'cArrowSurfacePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'cArrowSurfacePlotGroup'));
        obj.m.result.export('plotExporter2d').run();

        obj.m.result('cArrowLinePlotGroup').set('data',dset);
        obj.m.result('cArrowLinePlotGroup').set('view', 'standardView');
        obj.m.result('cArrowLinePlotGroup').feature('cArrowLinePlot').set('expr', {obj.cx_id{i} obj.cy_id{i}});
        obj.m.result('cArrowLinePlotGroup').feature('cArrowLinePlot').set('arrowlength', arrowLength);
        obj.m.result('cArrowLinePlotGroup').feature('cArrowLinePlot').set('scale', arrowScale);
        obj.m.result('cArrowLinePlotGroup').feature('cArrowLinePlot').set('descr', 'Gradient of c');
        obj.m.result('cArrowLinePlotGroup').feature('cArrowLinePlot').set('scaleactive', arrowScaleActive);
        obj.m.result.export('plotExporter2d').set('plotgroup', 'cArrowLinePlotGroup');
        obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'cArrowLinePlotGroup'));
%         obj.m.result('cArrowLinePlotGroup').feature('cArrowLinePlot').set('wireframecolor', 'none');
        obj.m.result.export('plotExporter2d').run();
    end
    
    %% current plot
    fileNameTemplate = [strrep(standard2dSubFolder,'\','\\'),'\\%s.png'];

    obj.m.result('iStreamlinePlotGroup').set('data',dset);
    obj.m.result('iStreamlinePlotGroup').set('view', 'standardView');
    obj.m.result('iStreamlinePlotGroup').feature('iStreamlinePlot').set('descr', 'current density');
    obj.m.result('iStreamlinePlotGroup').feature('iStreamlinePlot').set('posmethod', streamlinePosmeth);
    obj.m.result('iStreamlinePlotGroup').feature('iStreamlinePlot').set('udist', streamlineDist);
%     obj.m.result('iStreamlinePlotGroup').feature('iStreamlinePlot').set('mdensity', '15');
    obj.m.result('iStreamlinePlotGroup').feature('iStreamlinePlot').set('expr', {'ix' 'iy'});
    obj.m.result.export('plotExporter2d').set('plotgroup', 'iStreamlinePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'iStreamlinePlotGroup'));
    obj.m.result.export('plotExporter2d').run();

    obj.m.result('iArrowSurfacePlotGroup').set('data',dset);
    obj.m.result('iArrowSurfacePlotGroup').set('view', 'standardView');
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('xnumber', '20');
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('ynumber', '5');
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('expr', {'ix' 'iy'});
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('arrowlength', arrowLength);
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('scale', arrowScale);
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('descr', 'current density');
    obj.m.result('iArrowSurfacePlotGroup').feature('iArrowSurfacePlot').set('scaleactive', arrowScaleActive);
    obj.m.result.export('plotExporter2d').set('plotgroup', 'iArrowSurfacePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'iArrowSurfacePlotGroup'));
    obj.m.result.export('plotExporter2d').run();

    obj.m.result('iArrowLinePlotGroup').set('data',dset);
    obj.m.result('iArrowLinePlotGroup').set('view', 'standardView');
    obj.m.result('iArrowLinePlotGroup').feature('iArrowLinePlot').set('expr', {'ix' 'iy'});
    obj.m.result('iArrowLinePlotGroup').feature('iArrowLinePlot').set('arrowlength', arrowLength);
    obj.m.result('iArrowLinePlotGroup').feature('iArrowLinePlot').set('scale', arrowScale);
    obj.m.result('iArrowLinePlotGroup').feature('iArrowLinePlot').set('descr', 'Gradient of c');
    obj.m.result('iArrowLinePlotGroup').feature('iArrowLinePlot').set('scaleactive', arrowScaleActive);
    obj.m.result.export('plotExporter2d').set('plotgroup', 'iArrowLinePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'iArrowLinePlotGroup'));
    obj.m.result.export('plotExporter2d').run();
    
    %% conductivity plots
    obj.m.result('kappaContourPlotGroup').set('data',dset);
    obj.m.result('kappaContourPlotGroup').set('view', 'standardView');
    obj.m.result('kappaContourPlotGroup').feature('kappaContourPlot').label('kappa');
    obj.m.result('kappaContourPlotGroup').feature('kappaContourPlot').set('descr', 'Conductivity kappa');
    obj.m.result('kappaContourPlotGroup').feature('kappaContourPlot').set('expr', 'kappa_loc');
    obj.m.result.export('plotExporter2d').set('plotgroup', 'kappaContourPlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'kappaContourPlotGroup'));
    obj.m.result.export('plotExporter2d').run();
    
    obj.m.result('kappaSurfacePlotGroup').set('data',dset);
    obj.m.result('kappaSurfacePlotGroup').set('view', 'standardView');
    obj.m.result('kappaSurfacePlotGroup').feature('kappaSurfacePlot').label('phi');
    obj.m.result('kappaSurfacePlotGroup').feature('kappaSurfacePlot').set('descr', 'Dependent variable phi');
    obj.m.result('kappaSurfacePlotGroup').feature('kappaSurfacePlot').set('expr', 'kappa_loc');
    obj.m.result.export('plotExporter2d').set('plotgroup', 'kappaSurfacePlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'kappaSurfacePlotGroup'));
    obj.m.result.export('plotExporter2d').run();
    
    %% mesh plot
%     fileNameTemplate = [fullProjectPath,'\\%s.png'];

    obj.m.result('meshPlotGroup').set('data',dset);
    obj.m.result('meshPlotGroup').set('view', 'standardView');
    obj.m.result('meshPlotGroup').feature('meshPlot').set('wireframecolor', 'none');
    obj.m.result.export('plotExporter2d').set('plotgroup', 'meshPlotGroup');
    obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'meshPlotGroup'));
    obj.m.result.export('plotExporter2d').run();

end