function obj = plotStandard2d(obj,dset,plots)
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
    fileNameTemplate = [strrep(standard2dSubFolder,'\','\\'),'\\%s_%s.png'];
    
    [titles,expressions,desc] = extractPlotInfo(plots);
        
    nExpressions = numel(expressions);

    streamlinePosmeth = 'uniform'; % magnitude
    streamlineDist = 0.02;
    arrowScale = 0.01;
    arrowScaleActive = false;
%     arrowScaleActive = true;
%     arrowLength = 'normalized'; % 'logarithmic'
    arrowLength = 'proportional'; % 'logarithmic'
    
%% potential plots

    for i=1:nExpressions
        if ~iscell(expressions{i}) || numel(expressions{i}) == 1          
            obj.m.result('surfacePlotGroup').set('data',dset);
        %     obj.m.result('surfacePlotGroup').set('view', 'standardView');
            obj.m.result('surfacePlotGroup').feature('surfacePlot').set('descr', desc{i});
            obj.m.result('surfacePlotGroup').feature('surfacePlot').set('expr', expressions{i});
            obj.m.result('surfacePlotGroup').run;
            obj.m.result.export('plotExporter2d').set('plotgroup', 'surfacePlotGroup');
            obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,titles{i},'surfacePlotGroup'));
            obj.m.result.export('plotExporter2d').run();

            obj.m.result('contourPlotGroup').set('data',dset);
        %     obj.m.result('contourPlotGroup').set('view', 'standardView');
            obj.m.result('contourPlotGroup').feature('contourPlot').set('descr', desc{i});
            obj.m.result('contourPlotGroup').feature('contourPlot').set('expr', expressions{i});
            obj.m.result.export('plotExporter2d').set('plotgroup', 'contourPlotGroup');
            obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,titles{i},'contourPlotGroup'));
            obj.m.result.export('plotExporter2d').run();
        else % >= 2
%             curExp = expressions{i};
            obj.m.result('streamlinePlotGroup').set('data',dset);
        %     obj.m.result('streamlinePlotGroup').set('view', 'standardView');
            obj.m.result('streamlinePlotGroup').feature('streamlinePlot').set('descr', 'Gradient of phi');
            obj.m.result('streamlinePlotGroup').feature('streamlinePlot').set('posmethod', streamlinePosmeth);
            obj.m.result('streamlinePlotGroup').feature('streamlinePlot').set('udist', streamlineDist);
        %     obj.m.result('streamlinePlotGroup').feature('streamlinePlot').set('mdensity', '15');
            obj.m.result('streamlinePlotGroup').feature('streamlinePlot').set('expr', expressions{i});
            obj.m.result.export('plotExporter2d').set('plotgroup', 'streamlinePlotGroup');
            obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,titles{i},'streamlinePlotGroup'));
            obj.m.result.export('plotExporter2d').run();
   
            obj.m.result('arrowSurfacePlotGroup').set('data',dset);
%             obj.m.result('arrowSurfacePlotGroup').set('view', 'standardView');
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('xnumber', '20');
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('ynumber', '5');
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('expr', expressions{i});
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('arrowlength', arrowLength);
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('scale', arrowScale);
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('descr', 'Gradient of phi');
            obj.m.result('arrowSurfacePlotGroup').feature('arrowSurfacePlot').set('scaleactive', arrowScaleActive);
            obj.m.result.export('plotExporter2d').set('plotgroup', 'arrowSurfacePlotGroup');
            obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,titles{i},'arrowSurfacePlotGroup'));
            obj.m.result.export('plotExporter2d').run();
    
            obj.m.result('arrowLinePlotGroup').set('data',dset);
%             obj.m.result('arrowLinePlotGroup').set('view', 'standardView');
            obj.m.result('arrowLinePlotGroup').feature('arrowLinePlot').set('expr', expressions{i});
            obj.m.result('arrowLinePlotGroup').feature('arrowLinePlot').set('arrowlength', arrowLength);
            obj.m.result('arrowLinePlotGroup').feature('arrowLinePlot').set('scale', arrowScale);
            obj.m.result('arrowLinePlotGroup').feature('arrowLinePlot').set('descr', 'Gradient of phi');
            obj.m.result('arrowLinePlotGroup').feature('arrowLinePlot').set('scaleactive', arrowScaleActive);
        %     obj.m.result('obj.arrowLinePlotGroup').feature('arrowLinePlot').set('wireframecolor', 'none');
            obj.m.result.export('plotExporter2d').set('plotgroup', 'arrowLinePlotGroup');
            obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,titles{i},'arrowLinePlotGroup'));
            obj.m.result.export('plotExporter2d').run();
        end
    end
  
    %% mesh plot

%     obj.m.result('meshPlotGroup').set('data',dset);
%     obj.m.result('meshPlotGroup').set('view', 'standardView');
%     obj.m.result('meshPlotGroup').feature('meshPlot').set('wireframecolor', 'none');
%     obj.m.result.export('plotExporter2d').set('plotgroup', 'meshPlotGroup');
%     obj.m.result.export('plotExporter2d').set('pngfilename', sprintf(fileNameTemplate,'meshPlotGroup'));
%     obj.m.result.export('plotExporter2d').run();
% 
end