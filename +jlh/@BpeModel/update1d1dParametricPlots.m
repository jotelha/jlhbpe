function obj = update1d1dParametricPlots(obj,dset,plotDset)
    if( ~exist('plotDset','var'))
        plotDset = dset;
    end
%% setup exporter
%     obj.plotExporter1d.set('view', 'standardView');
    obj.plotExporter1d.set('printunit', 'mm');
    obj.plotExporter1d.set('webunit', 'px');
    obj.plotExporter1d.set('printheight', '90');
    obj.plotExporter1d.set('webheight', '600');
    obj.plotExporter1d.set('printwidth', '120');
    obj.plotExporter1d.set('webwidth', '800');
    obj.plotExporter1d.set('printlockratio', 'off'); %  off
    obj.plotExporter1d.set('weblockratio', 'off'); % off
    obj.plotExporter1d.set('printresolution', '300');
    obj.plotExporter1d.set('webresolution', '96');
    obj.plotExporter1d.set('size', 'manualprint'); % current, manualprint, manualweb
    obj.plotExporter1d.set('antialias', 'on');
    obj.plotExporter1d.set('zoomextents', 'on'); % off
    obj.plotExporter1d.set('title', 'off');
%     obj.plotExporter1d.set('legend', 'on');
    obj.plotExporter1d.set('legend', 'off');
    obj.plotExporter1d.set('logo', 'on');
    obj.plotExporter1d.set('options', 'on');
    obj.plotExporter1d.set('fontsize', '9');
    obj.plotExporter1d.set('customcolor', [1 1 1]);
    obj.plotExporter1d.set('background', 'color');
    obj.plotExporter1d.set('axes', 'on');
%     obj.plotExporter1d.set('grid', 'on'); % jlh, for 3d
    obj.plotExporter1d.set('qualitylevel', '92');
    obj.plotExporter1d.set('qualityactive', 'on'); % off
    obj.plotExporter1d.set('imagetype', 'png');
    
%     fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');
    fullProjectPath = [pwd(),'\',obj.projectPath];
%     fileNameTemplate = [fullProjectPath,'\\%s_%s.png'];
    % subfolder for dataset:
    dsetSubFolder = [fullProjectPath,'\',plotDset];
    if( ~exist( dsetSubFolder,'dir') )
        mkdir(dsetSubFolder);
    end
    
    %% setup datasets
%     datasets = {    'bpeSurfaceResults',...
%                     'entireSurfaceResults',...
%                     'zetaPlaneResults',...
%                     'bulkBoundaryResults',...
%                     'weResults',...
%                     'ceResults',...
%                     'centralCrossectionResults',...
%                     'centralDDLCrossectionResults',...
%                     'leftBpeEdgeCrossectionResults',...
%                     'rightBpeEdgeCrossectionResults',...
%                     'cathodeCrossectionResults',...
%                     'anodeCrossectionResults',...
%                     'halfDebyeLengthResults'};
    selections = { 'regionOfFirstDebyeLength', 'all'};
    nSelections = numel(selections);
    
    plotgroups = {  'standardPhiPlotGroup',...
                    'standardConcentrationsPlotGroup',...
                    'logConcentrationsPlotGroup' };
    nPlotgroups = numel(plotgroups);
    
    % lg = cell(1,nPlotgroups);
%     tl = {'phi / U_T', 'concentrations c / c_{ref}', 'logarithmic concentrations log( c / c_{ref} )'};
%     xl = {'x/L','x/L','x/L'};
%     yl = {'phi / U_T', 'c / c_{ref}', 'log( c / c_{ref} )'};
%     
%     f = figure();
%     set(f,'Position',600*[0 0 nCols nRows/obj.widthToHeight]);

%% get parameter info
    info = mphsolinfo(obj.m,'Dataset',dset,'NU','on');
    nParameters = size(info.solpar,1);
    nRuns = info.sizesolvals / nParameters;
    
    parameterNames = {info.solpar};
    parameterValues = reshape(info.solvals,nParameters,nRuns)';
  
%     obj.updateDatasets(dset);
    obj.updatePlotsDimensionless(plotDset);

%% loop over datsets, plot groups and parameters
    for j=1:nSelections 
        selectionSubFolder = [dsetSubFolder,'\',selections{j}];
        if( ~exist( selectionSubFolder,'dir') )
            mkdir(selectionSubFolder);
        end
        
        for k=1:nPlotgroups
            plotgroupSubFolder = [selectionSubFolder,'\',plotgroups{k}];
            if( ~exist( plotgroupSubFolder,'dir') )
                mkdir(plotgroupSubFolder);
            end
             for l = 1:nRuns
                parameterSuffix = '';
                for m = 1:nParameters
                    parameterSuffix = ['_',parameterNames{m},'_',num2str(parameterValues(l,m)),parameterSuffix];
                end
                
                pathTemplate = strrep(plotgroupSubFolder,'\','\\');
                fileNameTemplate = [pathTemplate,'\\%s_%s',parameterSuffix,'.png'];
                
                % for all plots in plot group
                plots = char(obj.m.result(plotgroups{k}).feature.tags);
                nPlots = size(plots,1);
                for m = 1:nPlots
                    plot = strtrim( plots(m,:) );
                    if strcmp(selections{j},'all')
                        obj.m.result(plotgroups{k}).feature(plot).selection.all;
                    else
                        obj.m.result(plotgroups{k}).feature(plot).selection.named(selections{j});
                    end
                end
                obj.m.result(plotgroups{k}).set('solnum',l);
                obj.plotExporter1d.set('plotgroup', plotgroups{k});
                fileName = sprintf(fileNameTemplate,plotgroups{k},selections{j});
                obj.plotExporter1d.set('pngfilename', fileName);
                fprintf('  Saving plot as "%s"...\n', fileName);
                obj.plotExporter1d.run();
             end
        end
    end    
end
              
                    