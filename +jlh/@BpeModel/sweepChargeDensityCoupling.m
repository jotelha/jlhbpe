%% make a new project
function obj = sweepChargeDensityCoupling(obj)     
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    
    close all;

    obj.setKrawiecCaseGeneralParameters;
    obj.setKrawiecCaseReactionParameters;
    E_m = obj.calcMixedPotential;
    obj.plotCurrents;
    
    % E_m = phi_bpe;
    obj.phi_bpe = E_m;
    obj.phi0 = E_m / UT;
    % plotsVisible = 'on';
    % bpePoissonEquationBC = 'Dirichlet'; % no Stern layer

    obj.prepareIdentifiers();
    obj.createModel();
    
%% sweep

    obj.updateModel(); % needs to be called to set parameters

    parameterName = 'chargeDensityRampFactor';
    parameterLegend = 'coupling ramp factor';
    parameterValues = 0.01:0.01:1;

    obj.addParametricSweepStudy(parameterName,parameterLegend,parameterValues);
    %updateParametricSweepStudy
% 
%     fileName = makeFileName('pecletStabilityMesh.png');
%     meshManually; % needs to be called after adding a new study
%     clear fileName;

    obj.stationaryStudy.run();

    %% process results
    % selects latest dataset:
        dset = getLatestDataset(m);

    % updateProbes
    updatePlotsDimensionless;

    probeSurface;
    fileName = makeFileName('chargeDensityCouplingSweepSurfaceQuantities.png');
    plotParametricSweep;
    clear fileName;

    % % *.png files should be removed from ./dat before replacing by new plots
    fileName = makeFileName('chargeDensityCouplingSweepOverview.png');
    plotPotentialAndConcentrationOverview;
    clear fileName;
 % stores plots to file
% plotMixedPotentialConvergence
% plotSurfaceFluxConvergence
% plotSurfaceFluxLogConvergence
