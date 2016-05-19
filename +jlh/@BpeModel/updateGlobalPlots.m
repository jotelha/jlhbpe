function updateGlobalPlots(obj,dset,par)
    import jlh.*
    plotNames = {'i_total', 'i','log_i'};
    integrationTemplate = 'integrateSurface(%s)';
    expressions = { 'integrateSurface(i_total)', ...
                    {'integrateSurface(i_cathodic)','integrateSurface(i_anodic)'},...
                    {'log(abs(integrateSurface(i_cathodic)))','log(abs(integrateSurface(i_anodic)))'}};
    
    % plots for single reactions
    for i = 1:obj.nReactions
        plotNames{end+1} = obj.i_id{i};
        expressions{end+1} = sprintf(integrationTemplate,obj.i_id{i});
    end
    for i = 1:obj.numberOfSpecies
        plotNames{end+1} = obj.N_id{i};
        expressions{end+1} = sprintf(integrationTemplate,obj.N_id{i});
        
        plotNames{end+1} = obj.ny_id{i};
        expressions{end+1} = sprintf(integrationTemplate,obj.ny_id{i});
    end
        
        
    nExpressions = numel(expressions);
    fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');
    %fileNameTemplate = [fullProjectPath,'\\%s_%s.png'];
    
    
                
    info = mphsolinfo(obj.m,'Dataset',dset,'NU','on');
    nParameters = size(info.solpar,1);
    nRuns = info.sizesolvals / nParameters;
    
    parameterNames = {info.solpar};
    parameterValues = reshape(info.solvals,nParameters,nRuns)';
    
    plotParameterPosition = find( strcmp(parameterNames,par) );
    
%     looplevelinput = cell(1,nParameters);
    looplevelinput = cellstr(repmat('manual',nParameters,1))';
    looplevelinput{plotParameterPosition} = 'all';
    
    singleParameterValues = cell(1,nParameters);
    parameterValueCount = zeros(1,nParameters);
    
    for i = 1:nParameters
        singleParameterValues{i} = unique(parameterValues(:,i));
        parameterValueCount(i) = numel(singleParameterValues{i});
    end
    
%     rowPos = 1;
    currentLooplevelSetting = cellstr(repmat('1',nParameters,1))';
    currentLooplevelSetting{plotParameterPosition} = ...
        strjoin( strtrim(cellstr(num2str( ...
            (1:parameterValueCount(plotParameterPosition))')))',',');

    cellPos = 1;
%     increment = 1;   
    rowPos = ones(1,nParameters);
    cfgCount = 1;
    while(cellPos > 0)
        if( cellPos ~= plotParameterPosition )
            currentLooplevelSetting{cellPos} = num2str(rowPos(cellPos));
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
                msg = prepTerm('parameter: value', 'parameter', 'value', parameterNames, currentLooplevelSetting);
                cellfun( @(s) fprintf('%s\n',s), msg);
                cfgCount = cfgCount + 1;
                
                parameterSuffix = '';
                for i = 1:nParameters
                    if i ~= plotParameterPosition
                        val = singleParameterValues{i};
                        parameterSuffix = ['_',parameterNames{i},'_',num2str(val(rowPos(i)-1)),parameterSuffix];
                    end
                end
                fileNameTemplate = [fullProjectPath,'\\%s',parameterSuffix,'.png'];

                
                for i = 1:nExpressions
                    % obj.globalPlotGroup;
                    obj.globalPlotGroup.set('data',dset);

                    %obj.globalPlotGroup.label('1D Plot Group 3');
                    obj.globalPlotGroup.set('xlabel', par);
                    obj.globalPlotGroup.set('xlabelactive', false);
                    obj.globalPlotGroup.set('looplevelinput', fliplr(looplevelinput));
                    %obj.globalPlot.set('descr', {'' ''});
                    %obj.globalPlot.set('unit', {'m^2*A/mol' ''});
                    obj.globalPlot.set('looplevel',fliplr(currentLooplevelSetting));
                    obj.globalPlot.set('expr', expressions{i});
                    
                    obj.plotExporter1d.set('plotgroup', 'globalPlotGroup');
                    fileName = sprintf(fileNameTemplate,plotNames{i});
                    obj.plotExporter1d.set('pngfilename', fileName);
                    fprintf('  Saving plot as "%s"...\n', fileName);
                    obj.plotExporter1d.run();
                end

            end
        else          
            rowPos(cellPos) = 1;
            cellPos = cellPos - 1;         
        end
        
        
    end
            
%             
%             
% 
%         end
%         
%         if rowPos < parameterValueCount(cellPos)
%                 rowPos = 
%         
%            
%         if(cellPos == plotParameterPosition)
%             cellPos = cellPos + increment;
%         end
%           
%         rowPos(cellPos) = rowPos(cellPos) + 1;
%         
%         currentLoopLevelSetting{cellPos} = num2str(rowPos(cellPos));
% 
%         if cellPos == nParameters
%             fprintf(currentLoopLevelSetting);
%         end
%         increment = 1;
%         if rowPos(cellPos) >= parameterValueCount(cellPos)
%             rowPos(cellPos) = 0;
%             increment = -1;
% %             cellPos = cellPos - 1;
% %             continue
%         end
%     end
%         
%         
%     for i = 1:nParameters
%         if i == plotParameterPosition
%             continue;
%         end
%         for j = 1:parameterValueCount(i)
%             
%     end
%             
    
%     m.m.sol('sol1').feature('s1').feature('pDef').getStringArray('plistarr')
% model.result('pg3').feature('glob1').set('data', 'dset1');
   
% obj.globalPlotGroup;
%     obj.globalPlotGroup.set('data',dset);
%     
%     obj.globalPlotGroup.label('1D Plot Group 3');
%     obj.globalPlotGroup.set('xlabel', 'phiSymmetryFactor');
%     obj.globalPlotGroup.set('xlabelactive', false);
%     obj.globalPlotGroup.set('looplevelinput', {'all' 'manual'});
%     %obj.globalPlot.set('descr', {'' ''});
%     %obj.globalPlot.set('unit', {'m^2*A/mol' ''});
% %     obj.globalPlot.set('looplevel', {'1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21' '2'});
%     obj.globalPlot.set('expr', {'integrateSurface(i_total)' ''});

    
% expression = [ 'phi0','phi','i_total','i_anodic','i_cathodic',c_id',N_id' ];
% surfaceC = cell(1,numberOfSpecies);
% surfaceN = cell(1,numberOfSpecies);
% % [surfacePhi0, surfacePhi, surfaceItotal, surfaceIanodic, surfaceIcathodic, surfaceC{:},surfaceN{:}] = mphevalpoint(m, expression, 'selection', 'surfaceNode','Dataset',dset);
% [surfacePhi0, surfacePhi, surfaceItotal, surfaceIanodic, surfaceIcathodic, surfaceC{:},surfaceN{:}] = mphevalpoint(m, expression, 'selection', 'surfaceNode','Dataset',dset,'outersolnum','all');
% potentialFieldPattern       = ' %12.4f ';
% currentFieldPattern         = ' %12.4e ';
% concentrationFieldPattern   = ' %12.4e ';
% fluxFieldPattern            = ' %12.4e ';
% 
% surfacePhi0_str = aprintf(potentialFieldPattern,surfacePhi0);
% surfacePhi_str = aprintf(potentialFieldPattern,surfacePhi);
% surfaceItotal_str = aprintf(currentFieldPattern,surfaceItotal);
% surfaceIanodic_str = aprintf(currentFieldPattern,surfaceIanodic);
% surfaceIcathodic_str = aprintf(currentFieldPattern,surfaceIcathodic);
% surfaceC_str = cellfun(@(c) aprintf(concentrationFieldPattern,c), surfaceC, 'UniformOutput', false );
% surfaceN_str = cellfun(@(c) aprintf(fluxFieldPattern,c), surfaceN, 'UniformOutput', false );
% 
% % table head:
% %   phi0  phi  i_total  i_anodic  i_cathodic  c1 c2 ... cN  N1 N2 ... NN
% nFields = 5 + 2*numberOfSpecies;
% tableHeadFieldPattern = ' %12.12s ';
% tableHeadFields = [{'phi0', 'phi', 'i_total', 'i_anodic', 'i_cathodic'}, c_id', N_id'];
% tableHeadRowPattern = repmat(tableHeadFieldPattern,1,nFields);
% tableHead = sprintf(tableHeadRowPattern,tableHeadFields{:});
% tableRows = strcat(surfacePhi0_str,surfacePhi_str,surfaceItotal_str,surfaceIanodic_str,surfaceIcathodic_str,surfaceC_str{:},surfaceN_str{:});
% 
% fprintf('%s\n',tableHead,tableRows{:});
% 
% if exist('fileName','var') && ischar(fileName)
%     f = fopen(fileName,'w','n','UTF-8');
%     fprintf(f,'%s\n',tableHead,tableRows{:});
%     fclose(f);
% end
% % fprintf('Surface probe:\n');
% % fprintf('  phi0    = % 12.4e\n', surfacePhi0);
% % fprintf('  phi     = % 12.4e\n', surfacePhi);
% % fprintf('  i_total = % 12.4e\n', surfaceItotal);
% % fprintf('  i_anod  = % 12.4e\n', surfaceIanodic);
% % fprintf('  i_cath  = % 12.4e\n', surfaceIcathodic);
% 
% % for i=1:numberOfSpecies
% % fprintf('  %-8.8s= % 12.4e\n', c_id{i},surfaceC{i});
% % end
% % 
% % for i=1:numberOfSpecies
% % fprintf('  %-8.8s= % 12.4e\n', N_id{i},surfaceN{i});
% % end

end