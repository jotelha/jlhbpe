function obj = sampleProject(obj)
    import jlh.hf.*
    close all;
        
    obj.setKrawiecCaseGeneralParameters();
    obj.setKrawiecCaseReactionParameters();
    % calcMixedPotential
    % E_m = phi_bpe;

    obj.model_tag = 'sample_model';
    % plotsVisible = 'on';
    % bpePoissonEquationBC = 'Dirichlet'; % no Stern layer

    obj.prepareIdentifiers();
    obj.createModel();
%     attachStudy
% 
%     modifyModel % needs to be called to set parameters
%     %attachSolution;
%     addSolverSequence;
%     meshManually;
% 
%     stationaryStudy.run();
% 
%     % selects lates dataset:
%     % datasetTags = m.result.dataset.tags;
%     % dset = char(datasetTags(end));
% 
%     % selects first dataset:
%     dset = 'dset1';
% 
%     updateProbes
%     updatePlotsDimensionless;
% 
%     model{1} = m; % model cell array is used by plotting commands and for saving

    % 
    % 
    % % *.png files should be removed from ./dat before replacing by new plots
%    plotPotentialAndConcentrationOverview % stores plots to file
    % plotMixedPotentialConvergence
    % plotSurfaceFluxConvergence
    % plotSurfaceFluxLogConvergence
end