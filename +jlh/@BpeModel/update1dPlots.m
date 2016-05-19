function obj = update1dPlots(obj)
%% setup exporter
%     obj.m.m.result.export('plotExporter1d') .set('view', 'standardView');
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
    obj.m.result.export('plotExporter1d').set('zoomextents', 'off'); % off
    obj.m.result.export('plotExporter1d').set('title', 'off');
%     obj.m.result.export('plotExporter1d') .set('legend', 'on');
    obj.m.result.export('plotExporter1d').set('legend', 'off');
    obj.m.result.export('plotExporter1d').set('logo', 'off');
    obj.m.result.export('plotExporter1d').set('options', 'on');
    obj.m.result.export('plotExporter1d').set('fontsize', '9');
    obj.m.result.export('plotExporter1d').set('customcolor', [1 1 1]);
    obj.m.result.export('plotExporter1d').set('background', 'color');
    obj.m.result.export('plotExporter1d').set('axes', 'on');
%     obj.m.result.export('plotExporter1d') .set('grid', 'on'); % jlh, for 3d
    obj.m.result.export('plotExporter1d').set('qualitylevel', '92');
    obj.m.result.export('plotExporter1d').set('qualityactive', 'on'); % off
    obj.m.result.export('plotExporter1d').set('imagetype', 'png');
    
%     fullProjectPath = strrep([pwd(),'\',obj.projectPath],'\','\\');
    fullProjectPath = [pwd(),'\',obj.projectPath];

%     dsetSubFolder = [fullProjectPath,'\',dset];
%     if( ~exist( dsetSubFolder,'dir') )
%         mkdir(dsetSubFolder);
%     end
    
    %% setup datasets
    datasets = {    'bpeSurfaceResults',...
                    'entireSurfaceResults',...
                    'zetaPlaneResults',...
                    'bulkBoundaryResults',...
                    'weResults',...
                    'ceResults',...
                    'centralCrossectionResults',...
                    'centralDDLCrossectionResults',...
                    'leftBpeEdgeCrossectionResults',...
                    'rightBpeEdgeCrossectionResults',...
                    'cathodeCrossectionResults',...
                    'anodeCrossectionResults',...
                    'halfDebyeLengthResults'};
    nDatasets = numel(datasets);
    
    plotgroups = {  'standardPhiPlotGroup',...
                    'standardConcentrationsPlotGroup',...
                    'logConcentrationsPlotGroup' };
    nPlotgroups = numel(plotgroups);
    
    % multipurpose plots
    titles = {'phix','phiy'};
    expressions = { {'phix'}, {'phiy'} };
    ylabel = { 'phi_x * L / U_T', 'phi_y * L / U_T' };
    labels = expressions;
    nExpressions = numel(expressions);
    
    % lg = cell(1,nPlotgroups);
%     tl = {'phi / U_T', 'concentrations c / c_{ref}', 'logarithmic concentrations log( c / c_{ref} )'};
%     xl = {'x/L','x/L','x/L'};
%     yl = {'phi / U_T', 'c / c_{ref}', 'log( c / c_{ref} )'};
%     
    nRows = nDatasets;
    nCols = nPlotgroups;
%     f = figure();
%     set(f,'Position',600*[0 0 nCols nRows/obj.widthToHeight]);

    
    for j=1:nRows
        dset = datasets{j};
        sol = char(obj.m.result.dataset(dset).getString('data'));
        
        solSubFolder = [fullProjectPath,'\',sol];
        if( ~exist( solSubFolder,'dir') )
            mkdir(solSubFolder);
        end
        
        standard1dSubFolder = [solSubFolder,'\standard1d'];
        if( ~exist( standard1dSubFolder,'dir') )
            mkdir(standard1dSubFolder);
        end
        
        dsetSubFolder = [standard1dSubFolder,'\',dset];
        if( ~exist( dsetSubFolder,'dir') )
            mkdir(dsetSubFolder);
        end
        
        switch dset
            case {  'bpeSurfaceResults',...
                    'entireSurfaceResults',...
                    'zetaPlaneResults',...
                    'bulkBoundaryResults' }
                obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').active(false);
%                 lg1 = {char(obj.standardPhiPlot.getString('legends'))};
                
            otherwise
                obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').active(true);
%                 lg1 = {char(obj.standardPhiPlot.getString('legends'));...
%                         char(obj.phiPBPlot.getString('legends'))};
        end
        obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').active(false);
        obj.updatePlotsDimensionless(dset);

%         for h=1:obj.numberOfSpecies
%             lg2{h} = char(obj.standardConcentrationsPlot{h}.getString('legends'));
%         end
%         for h=1:obj.numberOfSpecies
%             lg3{h} = char(obj.logConcentrationsPlot{h}.getString('legends'));
%         end
%         lg = {lg1,lg2,lg3};
          
        fileNameTemplate = [strrep(dsetSubFolder,'\','\\'),'\\%s_%s.png'];

        for k=1:nCols
            obj.m.result.export('plotExporter1d').set('plotgroup', plotgroups{k});
            fileName = sprintf(fileNameTemplate,plotgroups{k},datasets{j});
            obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
            fprintf('  Saving plot as "%s"...\n', fileName);
            obj.m.result(plotgroups{k}).run();
            obj.m.result.export('plotExporter1d').run();
%             i = (j-1)*nCols + k;
%             subplot(nRows,nCols,i);
%             mphplot(obj.m,plotgroups{k});
%             legend(lg{k});
%             title(tl{k});
%             xlabel(xl{k});
%             ylabel(yl{k});
        end
        
        %% flexible multipurpose plots
         obj.m.result.export('plotExporter1d').set('plotgroup', 'multiPurpose1dPlotGroup');

         for k=1:nExpressions        
            curExp = expressions{k};
            curLbl = labels{k};
            nCurExp = numel( curExp );
            
            obj.m.result('multiPurpose1dPlotGroup').set('ylabel', ylabel{k});
            obj.m.result('multiPurpose1dPlotGroup').set('xlabel', 'x / L');
            for l=1:obj.nMultiPurpose1dPlots
                if l <= nCurExp
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(true);
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).set('expr', curExp{l});
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).label(curLbl{l});
                else
                    obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{l}).active(false);
                end
            end
            obj.m.result('multiPurpose1dPlotGroup').run;

            fileName = sprintf(fileNameTemplate,titles{k},datasets{j});
            obj.m.result.export('plotExporter1d').set('pngfilename', fileName);
            fprintf('  Saving plot as "%s"...\n', fileName);
            obj.m.result.export('plotExporter1d').run();
        end
    end    
end