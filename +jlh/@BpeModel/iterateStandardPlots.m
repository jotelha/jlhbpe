function obj = iterateStandardPlots(obj)
    datasets = {    'bpeSurfaceResults',...
                    'entireSurfaceResults',...
                    'zetaPlaneResults',...
                    'bulkBoundaryResults',...
                    'weResults',...
                    'ceResults',...
                    'centralCrossectionResults',...
                    'centralDDLCrossectionResults' };
    nDatasets = numel(datasets);
    
    plotgroups = {  'standardPhiPlotGroup',...
                    'standardConcentrationsPlotGroup',...
                    'logConcentrationsPlotGroup' };
    nPlotgroups = numel(plotgroups);
    
    % lg = cell(1,nPlotgroups);
    tl = {'phi / U_T', 'concentrations c / c_{ref}', 'logarithmic concentrations log( c / c_{ref} )'};
    xl = {'x/L','x/L','x/L'};
    yl = {'phi / U_T', 'c / c_{ref}', 'log( c / c_{ref} )'};
    
    nRows = nDatasets;
    nCols = nPlotgroups;
    f = figure();
    set(f,'Position',600*[0 0 nCols nRows/obj.widthToHeight]);

    
    for j=1:nRows
        dset = datasets{j};
        obj.updatePlotsDimensionless(dset);
        switch dset
            case {  'bpeSurfaceResults',...
                    'entireSurfaceResults',...
                    'zetaPlaneResults',...
                    'bulkBoundaryResults' }
                obj.phiPBPlot.active(false);
                lg1 = {char(obj.standardPhiPlot.getString('legends'))};
                
            otherwise
                obj.phiPBPlot.active(true);
                lg1 = {char(obj.standardPhiPlot.getString('legends'));...
                        char(obj.phiPBPlot.getString('legends'))};
        end
        for h=1:obj.numberOfSpecies
            lg2{h} = char(obj.standardConcentrationsPlot{h}.getString('legends'));
        end
        for h=1:obj.numberOfSpecies
            lg3{h} = char(obj.logConcentrationsPlot{h}.getString('legends'));
        end
        lg = {lg1,lg2,lg3};
                
        for k=1:nCols
            i = (j-1)*nCols + k;
            subplot(nRows,nCols,i);
            mphplot(obj.m,plotgroups{k});
            legend(lg{k});
            title(tl{k});
%             xlabel(xl{k});
%             ylabel(yl{k});
        end
    end    
    obj.savePlot(f,'standardPlots2d.png',nRows,nCols);
end
              
                    