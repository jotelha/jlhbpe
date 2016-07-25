function plotParametric1d(obj,dset,dsetFolder,selections,plots,step)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.*
    import jlh.hf.*
    
    if ~exist('dsetFolder','var')
        dsetFolder = dset;
    end
    if ~iscell(selections)
        selections = {selections};
    end
    if ~exist('step','var')
        step = 1;
    end
    
    % extract plot info
    titles = plots.keys';
    
    values = plots.values';
    hasLabel = cellfun(@(c) iscell(c),values);
    ylabel = cell(numel(titles),1);
    ylabel(hasLabel) = cellfun(@(c) c{2}, values(hasLabel),'UniformOutput',false);
    values(hasLabel) = cellfun(@(c) c{1}, values(hasLabel),'UniformOutput',false);
    expressions = values;
    
    %% get parameter info
%     info = mphsolinfo(obj.m,'Dataset',dset,'NU','on');
%     nParameters = size(info.solpar,1);
%     nRuns = info.sizesolvals / nParameters;

%     parameterNames = {info.solpar};
%     parameterValues = reshape(info.solvals,nParameters,nRuns)';
    
    info = mphsolutioninfo(obj.m,'dataset',dset);
    sol_id = info.solutions{1};
    
    parameterNames = info.(sol_id).parameters;
    nParameters = size(parameterNames,1);
    
    nRuns = size(info.(sol_id).map,1);
    parameterValues = info.(sol_id).map(:,1:(end-2));

    
    
    %% setup exporter
    obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');
    obj.m.result('multiPurpose1dPlotGroup').set('data',dset);

    fullProjectPath = [pwd(),'\',obj.projectPath];


    % selections
    nSelections = numel(selections);
    
    labels = expressions;
    nExpressions = numel(expressions);
    
 
%     solution1dSubFolder = [fullProjectPath,'\solution1d'];
%     if( ~exist( solution1dSubFolder,'dir') )
%         mkdir(solution1dSubFolder);
%     end
        
%     dsetSubFolder = [solution1dSubFolder,'\',dsetFolder];
    dsetSubFolder = [fullProjectPath,'\',dsetFolder];
    if( ~exist( dsetSubFolder,'dir') )
        mkdir(dsetSubFolder);
    end
        
        
    for m=1:nSelections
        selectionSubFolder = [dsetSubFolder,'\',selections{m}];
        if( ~exist( selectionSubFolder,'dir') )
            mkdir(selectionSubFolder);
        end
%         fileNameTemplate = [strrep(selectionSubFolder,'\','\\'),'\\%s_%s.png'];
        selectionSubFolder = [dsetSubFolder,'\',selections{m}];

        for k=1:nExpressions        
            expressionSubFolder = [selectionSubFolder,'\',titles{k}];

            curExp = expressions{k};
            curLbl = labels{k};
            nCurExp = numel( curExp );

            obj.m.result('multiPurpose1dPlotGroup').set('ylabel', ylabel{k});
            obj.m.result('multiPurpose1dPlotGroup').set('xlabel', 'x / L');
            for l=1:obj.nMultiPurpose1dPlots
                if l <= nCurExp
                    % export raw data
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(true);
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('expr', curExp{l});
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).label(curLbl{l});
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('data', 'parent');
                    if strcmp('all',selections{m})
                        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).selection.all();
                    else
                        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).selection.named(selections{m});
                    end
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdata', 'expr'); 
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataexpr', 'x');
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataunit', '1');
                else
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(false);
                end
            end
            
            for l = 1:step:nRuns
                parameterSuffix = '';
                for j = 1:nParameters
                    parameterSuffix = ['_',parameterNames{j},'_',num2str(parameterValues(l,j)),parameterSuffix];
                end
                parameterSuffix = [ sprintf('_%04d',l), parameterSuffix ];
%                 pathTemplate = strrep(plotgroupSubFolder,'\','\\');
%                 fileNameTemplate = [pathTemplate,'\\%s_%s',parameterSuffix,'.png'];
                fileNameTemplate = [strrep(expressionSubFolder,'\','\\'),'\\%s_%s',parameterSuffix,'.png'];

                obj.m.result('multiPurpose1dPlotGroup').set('solnum',l);
                obj.m.result('multiPurpose1dPlotGroup').run;

                obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');
                fileName = sprintf(fileNameTemplate,titles{k},selections{m});
                obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
                fprintf('  Saving plot as "%s"...\n', fileName);
                obj.m.result.export('plotExporter1d').run();
                
                if k==1
                     % export data
                    csvFileNameTemplate = [strrep(dsetSubFolder,'\','\\'),'\\%s',parameterSuffix,'.csv'];
                    csvFileName = sprintf(csvFileNameTemplate,selections{m});

                    e = ['x',flattenCell(expressions)];
                    if strcmp('all',selections{m})
                            d = mpheval(obj.m,e,'dataset',dset,'dataonly','on','solnum',l);
                        else
                            d = mpheval(obj.m,e,'dataset',dset,'dataonly','on','solnum',l,'selection',selections{m});
                    end
                    v = arrayfun(@(f) (d.(sprintf('d%d',f)))', 1:numel(e), 'UniformOutput', false);
                    out = cell2mat(v);
                    fid=fopen(csvFileName,'w');           
                    fprintf(fid,'%s,',e{1:(end-1)});
                    fprintf(fid,'%s\n',e{end});
                    fclose(fid);
                    dlmwrite(csvFileName,out,'-append','delimiter',',','precision',8);
                end
            end
%             obj.m.result('multiPurpose1dPlotGroup').run;
% 
%             fileName = sprintf(fileNameTemplate,titles{k},dset);
%             obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
%             fprintf('  Saving plot as "%s"...\n', fileName);
%             obj.m.result.export('plotExporter1d').run();
        end
    end
end
              
                    