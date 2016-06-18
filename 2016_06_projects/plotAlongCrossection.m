function plotAlongCrossection(obj,dset,xexp,plots,dsetFolder)
    import jlh.*
    import jlh.hf.*
    
    if ~exist('dsetFolder','var')
        dsetFolder = dset;
    end
    if ~exist('xexp','var')
        xexp = 'x';
    end
    
    % extract plot info
    titles = plots.keys';
    
    values = parameters.values';
    hasLabel = cellfun(@(c) iscell(c),values);
    ylabel = cell(numel(parameter),1);
    ylabel(hasLabel) = cellfun(@(c) c{2}, values(hasLabel),'UniformOutput',false);
    values(hasLabel) = cellfun(@(c) c{1}, values(hasLabel),'UniformOutput',false);
    expressions = values;
    
%% setup exporter
    obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');
    obj.m.result('multiPurpose1dPlotGroup').set('data',dset);


    fullProjectPath = [pwd(),'\',obj.projectPath];

   
    
    labels = expressions;
    nExpressions = numel(expressions);
    
    % lg = cell(1,nPlotgroups);
%     tl = {'phi / U_T', 'concentrations c / c_{ref}', 'logarithmic concentrations log( c / c_{ref} )'};
%     xl = {'x/L','x/L','x/L'};
%     yl = {'phi / U_T', 'c / c_{ref}', 'log( c / c_{ref} )'};
%     
%     nRows = nDatasets;
%     nCols = nPlotgroups;
%     f = figure();
%     set(f,'Position',600*[0 0 nCols nRows/obj.widthToHeight]);

    
%     for j=1:nRows
%         dset = datasets{j};
        sol = char(obj.m.result.dataset(dset).getString('data'));
        
        solSubFolder = [fullProjectPath,'\',sol];
        if( ~exist( solSubFolder,'dir') )
            mkdir(solSubFolder);
        end
        
        enhanced1dSubFolder = [solSubFolder,'\enhanced1d'];
        if( ~exist( enhanced1dSubFolder,'dir') )
            mkdir(enhanced1dSubFolder);
        end
        
        dsetSubFolder = [enhanced1dSubFolder,'\',dsetFolder];
        if( ~exist( dsetSubFolder,'dir') )
            mkdir(dsetSubFolder);
        end
        
        
%         obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').active(false);
%         obj.updatePlotsDimensionless(dset);

        
        fileNameTemplate = [strrep(dsetSubFolder,'\','\\'),'\\%s_%s.png'];
        
        for k=1:nExpressions        
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
%                     obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdata', 'arc'); 
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataunit', '1');
                else
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(false);
                end
                obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdata', 'expr'); 
                obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataexpr', xexp);
            end
            obj.m.result('multiPurpose1dPlotGroup').run;

            fileName = sprintf(fileNameTemplate,titles{k},dset);
            obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
            fprintf('  Saving plot as "%s"...\n', fileName);
            obj.m.result.export('plotExporter1d').run();
            
            % export data
            fileName(end-2:end) = 'txt';
            e = ['x','y',cellfun(@(c) c{:},expressions,'UniformOutput',false)];
            v = cell(numel(e),1);
            [v{:}] = mphinterp(obj.m,e,'dataset',dset);
            out = cell2mat(v)';
%             save(fileName,'v','-ascii','-double','-tabs');
%             dlmwrite(fileName,e,'delimiter',';');
            fid=fopen(fileName,'w');           
            fprintf(fid,'%s;',e{1:(end-1)});
            fprintf(fid,'%s\n',e{end});
            fclose(fid);
            dlmwrite(fileName,out,'-append','delimiter',';','precision',8);
        end
end
              
                    