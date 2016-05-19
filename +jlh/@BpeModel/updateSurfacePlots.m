function obj = updateSurfacePlots(obj,dset,plotDset)
    if( ~exist('plotDset','var'))
        plotDset = dset;
    end
    fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');

    import jlh.*
%     plotNames = {'i_total', 'i_cathodic','i_anodic',...
%                     'log_i_cathodic', 'log_i_anodic'};
    plotNames = {'i_total', 'i', 'log_i'};
%     integrationTemplate = 'integrateSurface(%s)';
    expressions = { {'i_total'}, ...
                    {'i_cathodic', 'i_anodic'},...
                    {'log(abs(i_cathodic))','log(abs(i_anodic))'}};
%     labels =  { 'i_total', ...
%                 'i_cathodic and i_anodic',...
%                 'log(abs(i_cathodic)) and log(abs(i_anodic))'};
    labels = expressions;
                
    % plots for single reactions
    for i = 1:obj.nReactions
        plotNames{end+1} = obj.i_id{i};
        expressions{end+1} = {obj.i_id{i}};
        labels{end+1} = obj.i_id(i);
    end
    for i = 1:obj.numberOfSpecies
        plotNames{end+1} = obj.N_id{i};
        expressions{end+1} = {obj.N_id{i}};
        labels{end+1} = obj.N_id(i);
    end
    
    nExpressions = numel(expressions);
       
    dsetSubFolder = [fullProjectPath,'\',plotDset];
    if( ~exist( dsetSubFolder,'dir') )
        mkdir(dsetSubFolder);
    end
%     plotSubFolder = [dsetSubFolder,'\surface',];
%     if( ~exist( plotSubFolder,'dir') )
%         mkdir(plotSubFolder);
%     end
%% setup exporter
%     obj.m.result.export('plotExporter1d').set('view', 'standardView');
    obj.m.result.export('plotExporter1d').set('printunit', 'mm');
    obj.m.result.export('plotExporter1d').set('webunit', 'px');
    obj.m.result.export('plotExporter1d').set('printheight', '90');
    obj.m.result.export('plotExporter1d').set('webheight', '600');
    obj.m.result.export('plotExporter1d').set('printwidth', '120');
    obj.m.result.export('plotExporter1d').set('webwidth', '800');
    obj.m.result.export('plotExporter1d').set('printlockratio', 'off'); %  off
    obj.m.result.export('plotExporter1d').set('weblockratio', 'off'); % off
    obj.m.result.export('plotExporter1d').set('printresolution', '300');
    obj.m.result.export('plotExporter1d').set('webresolution', '96');
    obj.m.result.export('plotExporter1d').set('size', 'manualprint'); % current, manualprint, manualweb
    obj.m.result.export('plotExporter1d').set('antialias', 'on');
    obj.m.result.export('plotExporter1d').set('zoomextents', 'on'); % off
    obj.m.result.export('plotExporter1d').set('title', 'off');
%     obj.m.result.export('plotExporter1d').set('legend', 'on');
    obj.m.result.export('plotExporter1d').set('legend', 'off');
    obj.m.result.export('plotExporter1d').set('logo', 'on');
    obj.m.result.export('plotExporter1d').set('options', 'on');
    obj.m.result.export('plotExporter1d').set('fontsize', '9');
    obj.m.result.export('plotExporter1d').set('customcolor', [1 1 1]);
    obj.m.result.export('plotExporter1d').set('background', 'color');
    obj.m.result.export('plotExporter1d').set('axes', 'on');
%     obj.m.result.export('plotExporter1d').set('grid', 'on'); % jlh, for 3d
    obj.m.result.export('plotExporter1d').set('qualitylevel', '92');
    obj.m.result.export('plotExporter1d').set('qualityactive', 'on'); % off
    obj.m.result.export('plotExporter1d').set('imagetype', 'png');
    
%     fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');
%     fileNameTemplate = [fullProjectPath,'\\%s_%s.png'];

%     dsetSubFolder = [fullProjectPath,'\',plotDset];
%     if( ~exist( dsetSubFolder,'dir') )
%         mkdir(dsetSubFolder);
%     end
    
    %% setup selections
    selections = { 'geom_bpeSurface' };
%     selections = { 'bpeSurfaceResults' };
    nSelections = numel(selections);
    
%     plotgroups = {  'multiPurpose1dPlotGroup' };
%     nPlotgroups = numel(plotgroups);

    %% get parameter info
    info = mphsolinfo(obj.m,'Dataset',dset,'NU','on');
    nParameters = size(info.solpar,1);
    nRuns = info.sizesolvals / nParameters;
    
    parameterNames = {info.solpar};
    parameterValues = reshape(info.solvals,nParameters,nRuns)';
    
    obj.m.result('multiPurpose1dPlotGroup').set('data',plotDset);

    for j=1:nSelections
        selectionSubFolder = [dsetSubFolder,'\',selections{j}];
        if( ~exist( selectionSubFolder,'dir') )
            mkdir(selectionSubFolder);
        end
                     
        for k=1:nExpressions
            curExp = expressions{k};
            curLbl = labels{k};
            nCurExp = numel( curExp );
            for l=1:obj.nMultiPurpose1dPlots
                if l <= nCurExp
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(true);
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('expr', curExp{l});
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).label([selections{j},'_',curLbl{l}]);
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).selection.named(selections{j});
                else
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(false);
                end
            end
%             obj.m.result('multiPurpose1dPlotGroup').feature('multiPurpose1dPlot').label([selections{j},'_',labels{k}]);
%             obj.m.result('multiPurpose1dPlotGroup').feature('multiPurpose1dPlot').set('xdata', 'arc'); % get x data from an expression
%             obj.m.result('multiPurpose1dPlotGroup').feature('multiPurpose1dPlot').set('expr', expressions{k});
%             obj.m.result('multiPurpose1dPlotGroup').feature('multiPurpose1dPlot').set('expr', expressions{k});
%             obj.m.result('multiPurpose1dPlotGroup').feature('multiPurpose1dPlot').selection.named(selections{j});


            for l=1:nRuns
                parameterSuffix = '';
                for m = 1:nParameters
                    parameterSuffix = ['_',parameterNames{m},'_',num2str(parameterValues(l,m)),parameterSuffix];
                end
                pathTemplate = strrep(selectionSubFolder,'\','\\');
                fileNameTemplate = [pathTemplate,'\\%s_%s',parameterSuffix,'.png'];
                fileName = sprintf(fileNameTemplate,plotNames{k},selections{j});
                obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');
                obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
                fprintf('  Saving plot as "%s"...\n', fileName);
                obj.m.result('multiPurpose1dPlotGroup').run();
                obj.m.result.export('plotExporter1d').run();
            end
            
        end
    end    
end
              
                    