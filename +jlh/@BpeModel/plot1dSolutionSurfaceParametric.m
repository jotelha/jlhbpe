function obj = plot1dSolutionSurfaceParametric(obj,dset,par,dsetFolder)
if ~exist('dsetFolder','var')
    dsetFolder = dset;
end

import jlh.*
import jlh.hf.*
%% setup exporter
obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');
obj.m.result('multiPurpose1dPlotGroup').set('data',dset);


%     obj.m.m.result.export('plotExporter1d') .set('view', 'standardView');

%     fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');
fullProjectPath = [pwd(),'\',obj.projectPath];

%     dsetSubFolder = [fullProjectPath,'\',dset];
%     if( ~exist( dsetSubFolder,'dir') )
%         mkdir(dsetSubFolder);
%     end

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
%     nDatasets = numel(datasets);
    
%     plotgroups = {  'standardPhiPlotGroup',...
%                     'standardConcentrationsPlotGroup',...
%                     'logConcentrationsPlotGroup' };
%     nPlotgroups = numel(plotgroups);

    % selections
%     selections = {'geom1d_ddl1dCumulative_dom','geom1d_extendedDdl1dCumulative_dom','geom1d_space1dCumulative_dom','all'};
%     nSelections = numel(selections);
    
    % multipurpose plots
    titles = {  'phi','phix','_PHI','_PHIx',...
                'c','logc','cx',...
                '_C','_logC', '_Cx',...
                'nx','ix','kappa'  };
    logc = prepTerm('log(c_id)','c_id',obj.c_id);
    cx = prepTerm('c_idx','c_id',obj.c_id);
%     cy = prepTerm('c_idy','c_id',obj.c_id);
%     absc  = prepTerm('sqrt(cx^2+cy^2)','cx','cy',cx,cy);
    C = prepTerm('c_id*c_ref','c_id',obj.c_id);
    logC = prepTerm('log(c_id*c_ref)','c_id',obj.c_id);
    Cx = prepTerm('(c_idx*c_ref/L)','c_id',obj.c_id);
%     Cy = prepTerm('(c_idy*c_ref/L)','c_id',obj.c_id);
%     absC  = prepTerm('sqrt(Cx^2+Cy^2)','Cx','Cy',Cx,Cy);
    
%     absn = prepTerm('sqrt(nx_id^2+ny_id^2)','nx_id','ny_id',obj.nx_id,obj.ny_id);

    nakedExpressions = {     {'phi'}, {'phix'},...
                        {'phi*UT'}, {'phix*UT/L'},...
                        obj.c_id, logc, cx, ...
                        C,logC,Cx,...
                        obj.nx_id,{'ix'},{'kappa_loc'} };
                    
    expressions = cellfun(@(e) ...
        cellfun( @(f) sprintf('integrateSurface1d(%s)',f),...
            e, 'UniformOutput',false),...
            nakedExpressions, 'UniformOutput',false);
                    
%     expressions = obj.domainExpressions1d;
%     surfaceExpressions = obj.surfaceExpressions1d;
                        
                       
    ylabel = { 'phi / U_T', 'phi_x * L / U_T', ...
                'phi / V', 'phi_x / V m^{-1}',...
                'c / c_ref', 'log(c/c_ref)','c_x * L / c_ref',...
                'c / mol m^{-3}', 'log(c/mol m^{-3})', 'c_x / mol m^{-4}',...
                'nx','ix','kappa'};
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
%         sol = char(obj.m.result.dataset(dset).getString('data'));
        
%         solSubFolder = [fullProjectPath,'\',sol];
%         if( ~exist( solSubFolder,'dir') )
%             mkdir(solSubFolder);
%         end
        
        solution1dSubFolder = [fullProjectPath,'\solution1d'];
        if( ~exist( solution1dSubFolder,'dir') )
            mkdir(solution1dSubFolder);
        end
        
        dsetSubFolder = [solution1dSubFolder,'\',dsetFolder];
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
        
        
%         obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').active(false);
%         obj.updatePlotsDimensionless(dset);

        
%         csvFileNameTemplate = [strrep(dsetSubFolder,'\','\\'),'\\%s.csv'];
        
        singleParameterValues = cell(1,nParameters);
        parameterValueCount = zeros(1,nParameters);
        
        for i = 1:nParameters
            singleParameterValues{i} = unique(parameterValues(:,i));
            parameterValueCount(i) = numel(singleParameterValues{i});
        end

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
%                     pathTemplate = strrep(globalSubFolder,'\','\\');
                    fileNameTemplate = [strrep(globalSubFolder,'\','\\'),'\\%s',parameterSuffix,'.png'];
%                     fileNameTemplate = [strrep(globalSubFolder,'\','\\'),'\\%s_%s.png'];

                    solnum = find( all(preselectedSolnum,2) );
    %                 solnumStr = strjoin( strtrim(cellstr(num2str(solnum))'),',');
                    solnumStr = strtrim(cellstr(num2str(solnum))');
                    for i = 1:nExpressions
                        curExp = expressions{i};
%                         curLbl = labels{k};
                        nCurExp = numel( curExp );
                        
                        % obj.m.result('globalPlotGroup');
                        
                        obj.m.result('globalPlotGroup').set('data',dset);

                        %obj.m.result('globalPlotGroup').label('1D Plot Group 3');
                        obj.m.result('globalPlotGroup').set('xlabel', par);
                        obj.m.result('globalPlotGroup').set('ylabel', ylabel{i});

                        obj.m.result('globalPlotGroup').set('xlabelactive', false);
                        obj.m.result('globalPlotGroup').set('innerinput', 'manual');
                        obj.m.result('globalPlotGroup').set('solnum', solnumStr);
    %                     obj.m.result('globalPlotGroup').set('looplevelinput', fliplr(looplevelinput));
                        %obj.globalPlot.set('descr', {'' ''});
                        %obj.globalPlot.set('unit', {'m^2*A/mol' ''});
    %                     obj.globalPlot.set('looplevel',fliplr(currentLooplevelSetting));
%                         for j=1:nCurExp
                        for l=1:obj.nMultiPurpose1dPlots
                            if l <= nCurExp
                               obj.m.result('globalPlotGroup').feature('globalPlot').setIndex('expr', curExp{l},l-1);
                            else
                               obj.m.result('globalPlotGroup').feature('globalPlot').setIndex('expr', '',l-1);
                            end
                        end
                        obj.m.result('globalPlotGroup').run;

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
            
%         for m=1:nSelections
%             selectionSubFolder = [dsetSubFolder,'\',selections{m}];
%             if( ~exist( selectionSubFolder,'dir') )
%                 mkdir(selectionSubFolder);
%             end
%             fileNameTemplate = [strrep(selectionSubFolder,'\','\\'),'\\%s_%s.png'];
% 
%             
%             for k=1:nExpressions        
%                 curExp = expressions{k};
%                 curLbl = labels{k};
%                 nCurExp = numel( curExp );
% 
%                 obj.m.result('multiPurpose1dPlotGroup').set('ylabel', ylabel{k});
%                 obj.m.result('multiPurpose1dPlotGroup').set('xlabel', 'x / L');
%                 for l=1:obj.nMultiPurpose1dPlots
%                     if l <= nCurExp
%                         % export raw data
%                         obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(true);
%                         obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('expr', curExp{l});
%                         obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).label(curLbl{l});
%                         obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('data', 'parent');
%                         if strcmp('all',selections{m})
%                             obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).selection.all();
%                         else
%                             obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).selection.named(selections{m});
%                         end
%                         obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdata', 'arc'); 
%     %                     obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataexpr', 'x');
%     %                     obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('xdataunit', '1');
%                     else
%                         obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(false);
%                     end
%                 end
%                 obj.m.result('multiPurpose1dPlotGroup').run;
% 
%                 fileName = sprintf(fileNameTemplate,titles{k},dset);
%                 obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
%                 fprintf('  Saving plot as "%s"...\n', fileName);
%                 obj.m.result.export('plotExporter1d').run();
%             end
            
%             % export data
%             csvFileName = sprintf(csvFileNameTemplate,selections{m});
% 
% %             e = ['x',cellfun(@(c) {c{:}},expressions,'UniformOutput',false)];
%             e = ['x',flattenCell(expressions)];
% %             v = cell(numel(e),1);
%             if strcmp('all',selections{m})
%                     d = mpheval(obj.m,e,'dataset',dset,'dataonly','on');
%                 else
%                     d = mpheval(obj.m,e,'dataset',dset,'dataonly','on','selection',selections{m});
%             end
%             v = arrayfun(@(f) (d.(sprintf('d%d',f)))', 1:numel(e), 'UniformOutput', false);
% %             [v{:}] = mphinterp(obj.m,e,'coord',0,'dataset',dset);
%             out = cell2mat(v);
% %             save(fileName,'v','-ascii','-double','-tabs');
% %             dlmwrite(fileName,e,'delimiter',';');
%             fid=fopen(csvFileName,'w');           
%             fprintf(fid,'%s,',e{1:(end-1)});
%             fprintf(fid,'%s\n',e{end});
%             fclose(fid);
%             dlmwrite(csvFileName,out,'-append','delimiter',',','precision',8);
%         end
end
              
                    