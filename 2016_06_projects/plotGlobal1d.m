function plotGlobal1d(obj,dset,par,plots,plotDset)
    if( ~exist('plotDset','var'))
        plotDset = dset;
    end
    import jlh.*
%     plotNames = {'i_total', 'i','log_i'};
%     integrationTemplate = 'integrateSurface(%s)';
%     expressions = { 'integrateSurface(i_total)', ...
%                     {'integrateSurface(i_cathodic)','integrateSurface(i_anodic)'},...
%                     {'log(abs(integrateSurface(i_cathodic)))','log(abs(integrateSurface(i_anodic)))'}};
                 % extract plot info
    titles = plots.keys';
    values = plots.values';
    hasLabel = cellfun(@(c) iscell(c),values);
    ylabel = cell(numel(titles),1);
    ylabel(hasLabel) = cellfun(@(c) c{2}, values(hasLabel),'UniformOutput',false);
    values(hasLabel) = cellfun(@(c) c{1}, values(hasLabel),'UniformOutput',false);
    expressions = values;
%    %% settings for exporter:
%     obj.m.result.export('plotExporter1d').set('printunit', 'mm');
%     obj.m.result.export('plotExporter1d').set('webunit', 'px');
%     obj.m.result.export('plotExporter1d').set('printheight', '90');
%     obj.m.result.export('plotExporter1d').set('webheight', '600');
%     obj.m.result.export('plotExporter1d').set('printwidth', '120');
%     obj.m.result.export('plotExporter1d').set('webwidth', '800');
%     obj.m.result.export('plotExporter1d').set('printlockratio', 'off'); %  off
%     obj.m.result.export('plotExporter1d').set('weblockratio', 'off'); % off
%     obj.m.result.export('plotExporter1d').set('printresolution', '300');
%     obj.m.result.export('plotExporter1d').set('webresolution', '96');
%     obj.m.result.export('plotExporter1d').set('size', 'manualprint'); % current, manualprint, manualweb
%     obj.m.result.export('plotExporter1d').set('antialias', 'on');
%     obj.m.result.export('plotExporter1d').set('zoomextents', 'on'); % off
%     obj.m.result.export('plotExporter1d').set('title', 'off');
%     %     obj.m.result.export('plotExporter1d').set('legend', 'on');
%     obj.m.result.export('plotExporter1d').set('legend', 'off');
%     obj.m.result.export('plotExporter1d').set('logo', 'on');
%     obj.m.result.export('plotExporter1d').set('options', 'on');
%     obj.m.result.export('plotExporter1d').set('fontsize', '9');
%     obj.m.result.export('plotExporter1d').set('customcolor', [1 1 1]);
%     obj.m.result.export('plotExporter1d').set('background', 'color');
%     obj.m.result.export('plotExporter1d').set('axes', 'on');
%     %     obj.m.result.export('plotExporter1d').set('grid', 'on'); % jlh, for 3d
%     obj.m.result.export('plotExporter1d').set('qualitylevel', '92');
%     obj.m.result.export('plotExporter1d').set('qualityactive', 'on'); % off
%     obj.m.result.export('plotExporter1d').set('imagetype', 'png');
    
    % plots for single reactions
%     for i = 1:obj.nReactions
%         plotNames{end+1} = obj.i_id{i};
%         expressions{end+1} = sprintf(integrationTemplate,obj.i_id{i});
%     end
%     for i = 1:obj.numberOfSpecies
%         plotNames{end+1} = obj.N_id{i};
%         expressions{end+1} = sprintf(integrationTemplate,obj.N_id{i});
%         
% %         plotNames{end+1} = obj.ny_id{i};
% %         expressions{end+1} = sprintf(integrationTemplate,obj.ny_id{i});
%     end
        
        
    nExpressions = numel(expressions);
%     fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');
    %fileNameTemplate = [fullProjectPath,'\\%s_%s.png'];
    fullProjectPath = [pwd(),'\',obj.projectPath];
%     fileNameTemplate = [fullProjectPath,'\\%s_%s.png'];
    % subfolder for dataset:
    dsetSubFolder = [fullProjectPath,'\',plotDset];
    if( ~exist( dsetSubFolder,'dir') )
        mkdir(dsetSubFolder);
    end
    globalSubFolder = [dsetSubFolder,'\global'];
    if( ~exist( globalSubFolder,'dir') )
        mkdir(globalSubFolder);
    end
                   
    info = mphsolinfo(obj.m,'Dataset',dset,'NU','on');
    nParameters = size(info.solpar,1);
    nRuns = info.sizesolvals / nParameters;
    
    parameterNames = info.solpar; % with or without ' ?
    parameterValues = reshape(info.solvals,nParameters,nRuns)';
    
    plotParameterPosition = find( strcmp(parameterNames,par) );
    
    if isempty(plotParameterPosition)
    	fprintf('No such parameter!\n');
        return;
    end

        
    
%     looplevelinput = cell(1,nParameters);
%     looplevelinput = cellstr(repmat('manual',nParameters,1))';
%     looplevelinput{plotParameterPosition} = 'all';
%     looplevelinput{ strcmp(parameterNames,par) } = 'all';
    
    singleParameterValues = cell(1,nParameters);
    parameterValueCount = zeros(1,nParameters);
    
    for i = 1:nParameters
        singleParameterValues{i} = unique(parameterValues(:,i));
        parameterValueCount(i) = numel(singleParameterValues{i});
    end
    
%     rowPos = 1;
%     currentLooplevelSetting = cellstr(repmat('1',nParameters,1))';
%     currentLooplevelSetting{plotParameterPosition} = ...
%         strjoin( strtrim(cellstr(num2str( ...
%             (1:parameterValueCount(plotParameterPosition))')))',',');

    cellPos = 1;
%     increment = 1;   
    rowPos = ones(1,nParameters);
    cfgCount = 1;
    while(cellPos > 0)
        if( cellPos ~= plotParameterPosition )
%             currentLooplevelSetting{cellPos} = num2str(rowPos(cellPos));
        end
    
        if rowPos(cellPos) <= parameterValueCount(cellPos)
            rowPos(cellPos) = rowPos(cellPos) + 1;
            if( cellPos == plotParameterPosition )
                    rowPos(cellPos) = parameterValueCount(cellPos)+1;
            end
                
            if cellPos < nParameters
                cellPos = cellPos + 1;   
            else
                fprintf('Config %d...\n',cfgCount);
%                 msg = prepTerm('parameter: value', 'parameter', 'value', parameterNames, currentLooplevelSetting);
%                 cellfun( @(s) fprintf('%s\n',s), msg);
                cfgCount = cfgCount + 1;
                
                parameterSuffix = '';
                preselectedSolnum = ones(nRuns,nParameters);
                for i = 1:nParameters
                    if i ~= plotParameterPosition
                        val = singleParameterValues{i};
                        parameterSuffix = ['_',parameterNames{i},'_',num2str(val(rowPos(i)-1)),parameterSuffix];
                        preselectedSolnum(:,i) = ( parameterValues(:,i) == val(rowPos(i)-1) );
                    end
                end
                pathTemplate = strrep(globalSubFolder,'\','\\');
                fileNameTemplate = [pathTemplate,'\\%s',parameterSuffix,'.png'];

                solnum = find( all(preselectedSolnum,2) );
                solnumStr = strtrim(cellstr(num2str(solnum))');
                
                for i = 1:nExpressions
                    % obj.m.result('globalPlotGroup');
                    obj.m.result('globalPlotGroup').set('ylabel', ylabel{i});

                    obj.m.result('globalPlotGroup').set('data',plotDset);

                    obj.m.result('globalPlotGroup').set('xlabel', par);
                    obj.m.result('globalPlotGroup').set('xlabelactive', false);
                    obj.m.result('globalPlotGroup').set('innerinput', 'manual');
                    obj.m.result('globalPlotGroup').set('solnum', solnumStr);
                    obj.m.result('globalPlotGroup').feature('globalPlot').set('expr', expressions{i});
                    
                    obj.m.result.export('plotExporter1d').set('plotgroup', 'globalPlotGroup');
                    fileName = sprintf(fileNameTemplate,titles{i});
                    obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
                    fprintf('  Saving plot as "%s"...\n', fileName);
                    obj.m.result.export('plotExporter1d').run();
                end

            end
        else          
            rowPos(cellPos) = 1;
            cellPos = cellPos - 1;         
        end
        
        
    end
end