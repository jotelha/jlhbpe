function plotBoundarySelection(obj,dset,selection,xexp,plots,plotDset)
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
	[titles,expressions,ylabels] = extractPlotInfo(plots);
    plotlabels = expressions;
        
    nExpressions = numel(expressions);
    fullProjectPath = [pwd(),'\',obj.projectPath];

    % subfolder for dataset:
%     solution1dSubFolder = [fullProjectPath,'\solution1d'];
%     if( ~exist( solution1dSubFolder,'dir') )
%         mkdir(solution1dSubFolder);
%     end
%     dsetSubFolder = [solution1dSubFolder,'\',plotDset];
    dsetSubFolder = [fullProjectPath,'\',plotDset];
    if( ~exist( dsetSubFolder,'dir') )
        mkdir(dsetSubFolder);
    end
    surfaceSubFolder = [dsetSubFolder,'\surface'];
    if( ~exist( surfaceSubFolder,'dir') )
        mkdir(surfaceSubFolder);
    end
         
    fileNameTemplate = [strrep(surfaceSubFolder,'\','\\'),'\\%s.png'];
    
    obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');
    obj.m.result('multiPurpose1dPlotGroup').set('data',dset);

    for k=1:nExpressions        
            curExp = expressions{k};
            curLbl = plotlabels{k};
            nCurExp = numel( curExp );
            
            obj.m.result('multiPurpose1dPlotGroup').set('ylabel', ylabels{k});
            obj.m.result('multiPurpose1dPlotGroup').set('xlabel', 'x / L');
            for l=1:obj.nMultiPurpose1dPlots
                if l <= nCurExp
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(true);
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).selection.named(selection);
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('expr', curExp{l});
%                   obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).label(curLbl{l});
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('data', 'parent');
%                     obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdata', 'arc'); 
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataunit', '1');
                else
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(false);
                end
                obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdata', 'expr'); 
                obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataexpr', xexp);
            end
            obj.m.result('multiPurpose1dPlotGroup').run;

            fileName = sprintf(fileNameTemplate,titles{k});
            obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
            fprintf('  Saving plot as "%s"...\n', fileName);
            obj.m.result.export('plotExporter1d').run();
            
            % export data
%             fileName(end-2:end) = 'txt';
%             e = ['x','y',cellfun(@(c) c{:},expressions,'UniformOutput',false)];
%             v = cell(numel(e),1);
%             [v{:}] = mphinterp(obj.m,e,'dataset',dset);
%             out = cell2mat(v)';
% %             save(fileName,'v','-ascii','-double','-tabs');
% %             dlmwrite(fileName,e,'delimiter',';');
%             fid=fopen(fileName,'w');           
%             fprintf(fid,'%s;',e{1:(end-1)});
%             fprintf(fid,'%s\n',e{end});
%             fclose(fid);
%             dlmwrite(fileName,out,'-append','delimiter',';','precision',8);
    end
end