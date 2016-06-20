function plotGlobal1d(obj,dset,par,plots,plotDset)
    if( ~exist('plotDset','var'))
        plotDset = dset;
    end
    import jlh.*

	% extract plot info
    % titles = plots.keys';
    % values = plots.values';
    % hasLabel = cellfun(@(c) iscell(c),values);
    % ylabel = cell(numel(titles),1);
    % ylabel(hasLabel) = cellfun(@(c) c{2}, values(hasLabel),'UniformOutput',false);
    % values(hasLabel) = cellfun(@(c) c{1}, values(hasLabel),'UniformOutput',false);
    % expressions = values;        
	[titles,expressions,ylabel] = extractPlotInfo(plots);
        
    nExpressions = numel(expressions);
    fullProjectPath = [pwd(),'\',obj.projectPath];

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

    singleParameterValues = cell(1,nParameters);
    parameterValueCount = zeros(1,nParameters);
    
    for i = 1:nParameters
        singleParameterValues{i} = unique(parameterValues(:,i));
        parameterValueCount(i) = numel(singleParameterValues{i});
    end
    
    cellPos = 1;
    rowPos = ones(1,nParameters);
    cfgCount = 1;
    while(cellPos > 0)
        if( cellPos ~= plotParameterPosition )
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