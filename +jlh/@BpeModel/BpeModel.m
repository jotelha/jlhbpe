classdef BpeModel < handle
    properties (Constant)
        F = jlh.Constants.FaradayConstant;
    end
    properties
        projectName
        projectPath
        model_tag
        reactions
        nReactions
        
        % potentials and currents
        E_m % mixed potential
       
        i
        i_anodic
        i_cathodic
        i_tot
        i_anodic_tot
        i_cathodic_tot
        
        % mesh properties
        vFirstDebyeLength
        nvFirstDebyeLength
        intFirstDebyeLength 
     
        vRemaining
        nvRemaining
        intRemaining

        vExtendedDdl
        nvExtendedDdl
        intExtendedDdl
        
        distributionFirstDebyeLength
        distributionFirstDebyeLengthStr
        distributionRemaining
        distributionRemainingStr
        distributionExtendedDdl
        distributionExtendedDdlStr

%% parameters
        transformation = 'none'; % apply no transormation to PNP system
        bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
        surfaceOverpotentialExpression = 'ZetaPotential';
        
        explicitElectrodeGeometry = true;
        
        bpeFluxOn = false; % switch species flux due to reactions on or off'
        ramp = false;
        rampSurfaceFlux = false;
        rampSurfacePotential = false;
        
        hMaxFactor = 1e-2; % size of smallest element compared to size allowed by Peclet criterion

        widthToHeight = 5/3; % width : height relation for single plots
        plotsVisible = 'on'; % 'off' stores plots, but does not show them as pop ups

        %
        % set up system parameters
        %

        % environment and physical constants
        epsilon = 0.01; % relation Debye length : simulation domain depth, lambdaD / L
        delta = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
%         gamma = 0.001; % relation Debye length : simulation domain width, lambdaD / W0
        beta  = 0.0009;
        
        % l = 1;
%         w % = W/L = epsilon*W/lambda_D
        w_bpe
        w_insulatorLeft
        w_insulatorRight
        w_cathode % 'WE'
        w_anode % 'CE'
        w_bulkLeft
        w_bulkRight
        
        wMeshRatio = 100; % ratio of w_mesh to l (l=1)
        elementSizeAtSurface;
        extendedDdlFactor = 3;
        
        epsilon_r = jlh.Constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

        T = 298.15; % K, equivalent to 25 deg C

        numberOfSpecies = 2;
        
        % bulk concentrations
        pH          % = 7;
        cH3Op_bulk  % = 10^(-pH)*10^3; %mol/m^3
        cOHm_bulk   % = cH3Op_bulk;
        cFe2p_bulk  % = 0;

        % solute species properties
        c_bulk  % = [cH3Op_bulk,cOHm_bulk]; % bulk concentrations
        z       % = [1,-1]; % charge number of solute species
        D       % = [jlh.Constants.DiffusivityH3Op,1e-9]; % diffusivity of solute species

        PHI_bpe     % = 0.1; % [phi] = V
        PHI_bulk    % = 0;
        
        deltaPHI                % potential difference from we to ce
        phiSymmetryFactor = 0.5 % symmetry of potential difference around bpe mixed potential

        % pure iron:
        % atomicWeightFe = 55.845 * 1.660539040e-27; %u -> kg, https://en.wikipedia.org/wiki/Iron, https://en.wikipedia.org/wiki/Atomic_mass_unit
        % densityFe = 6.98 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t., https://en.wikipedia.org/wiki/Iron
        % numberConcentrationFe = densityFe / atomicWeightFe;
        % molarConcentrationFe = numberConcentrationFe / jlh.Constants.AvogadroConstant;

        % concentrations (initial or constant)
        % cFe = molarConcentrationFe; % concentration of iron in solid form
        cFe = 1; %absorbed into rate constant
        % cFe2p = 0; % assume constant zero concentration of dissolved ion Fe2+ everywhere
        cH2O = 55560; % mol/m^3
        % saturation concentration of O2 in fresh water at 25 deg C and 1 bar
        cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
        %cO2 = 0.29  % 0.29mM = 0.29e-3 mol/l = 0.29 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html
        cH2 = 0; % assume no hydrogen in solution
        % cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html
        
        cH3Op   % = cH3Op_bulk;
        cOHm    % = cOHm_bulk;
        cFe2p   % = cFe2p_bulk;
        
        customParameters = struct;
                 
%% identifiers
           
        comp_id = 'comp1';
        
        % dependent on number of species
        speciesNames
        D_id
        z_id
        c_id
        cx_id
        cy_id
        C_id
        c_bulk_id 
        c_0_id
        c_pb_id
        cInterpolation_id
        
        lambdaCBulk_id
        surfaceFluxBC_id
        bulkFluxBC_id
        bulkConcentrationBC_id
        
        nx_id % dimensionless flux
        ny_id
        Nx_id % dimensional flux
        Ny_id
        
        ix_id % dimensionless current
        iy_id
        Ix_id % dimensional current
        Iy_id
        
        cSurfaceProbe_id
        NSurfaceProbe_id
        iSurfaceProbe_id
        
        standardConcentrationsPlot_id
        logConcentrationsPlot_id
        logConcentrationsAtSurfacePlot_id
        logPBConcentrationsAtSurfacePlot_id
        
        nu_id_pattern

        % dependent on number of surface reactions:
        reactionNames
        k0_id
        E0_id
        beta_id
        n_id
        i_id
        i_dimless_id
        nu_id
        N_id
        N_dimless_id
        
        multiPurpose1dPlot_id
    end
    properties (Dependent)
        RT; 
        UT;
        
        L
        w
        W
        Wbpe
        WinsulatorLeft
        WinsulatorRight
        Wcathode % 'WE'
        Wanode % 'CE'
        WbulkLeft
        WbulkRight
        
        w_mesh % size of one mesh chop
        nMeshChops
        nSegmentsPerChop
               
        PHI0 % dimensional
        phi0 % dimensionless, initial value
        phi_bpe % dimensionless, BPE potential
        
        delta_phi % dimensionless
        
        PHI_we % dimensional
        PHI_ce
        phi_we % dimensionless
        phi_ce
        
        ionicStrength
        lambdaD
        
        c_ref
        D_ref
        
        kappa % ionic conductivity after Newman 2012
        
        nMultiPurpose1dPlots
        
        % expressions for evaluation
        surfaceExpressions1d
        domainExpressions1d
        
        surfaceExpressions2d
        domainExpressions2d
    end
    properties (Transient)
%% COMSOL objects
        m % comsol model
       
        % analytical functions
        phiInterpolation

        % the analytical solution of Poisson-Boltzmann problem
        phi_pb
        phi_pbx
        c_pb
        cInterpolation % aslo for initial value interpolation
        
        % geometry
        geom
        gapLeft
        gapRight
        space   % remeining space above ddl
        spaceArray
        ddl     % the thin layer of width epsilon above surface
        leftEdgeOfBpe
        rightEdgeOfBpe
        
        leftEdgeOfCathode
        rightEdgeOfCathode
        leftEdgeOfAnode
        rightEdgeOfAnode
        
        % selections
%         surfaceNode 
%         zetaNode
%         bulkNode
        leftBoundaryOfSurface
        rightBoundaryOfSurface
        leftBoundaryOfZetaPlane
        rightBoundaryOfZetaPlane
        
        allBoundaries
        bpeSurface
        bulkBoundary
        zetaPlane
        lateralBoundaryAtFirstDebyeLength
        lateralBoundaryAtRemainingDomain
        lateralBoundary
        upperBoundary
        electrodes
        workingElectrode
        counterElectrode
        insulator
        insulatorAdjacentToBpe
        entireSurface
        reactingSurface
        
        regionOfFirstDebyeLength
        regionRemaining
        
        % local definitions
        projectReactionPlaneToSurface
        % projectSurfaceToBulkExit
        integrateSurface
        integrateBulkExit
        integrateWE
        integrateCE
        integrateDomain
        maximumOnDomain

        % variables
        domainVariables
        surfaceVariables
        bulkVariables
        upperBoundaryVariables
        electrodeVariables
        weVariables
        ceVariables

        % physics
        NernstPlanckEquation
        FluxAtSurface
        FluxAtElectrodes
        FluxAtBulkBoundary
        BulkConcentrationsBC
        
        PoissonEquation
        bpeChargeDensityBC
        bpePotentialBC
        bulkPotentialBC
        electrodePotentialBC
        bulkBoundaryPotentialGradientBC
        entireSurfacePotentialGradientBC
        insulatorSurfacePotentialGradientBC
        
        zeroSurfaceCurrent
        zeroSurfaceCurrentEquation
        
        electrodeCurrent
        electrodeCurrentEquation

        % mesh
        standardMesh
        meshingEdge % for 1d model
        bpeSurfaceMeshingEdge
        bpeSurfaceMeshingProperties
        
        zetaPlaneMeshingEdge
        zetaPlaneMeshingProperties
        
%         lateralBoundaryMeshingEdgeAtFirstDebyeLength
        lateralBoundaryMeshingPropertiesAtFirstDebyeLength
        lateralBoundaryMeshingPropertiesAtRemainingDomain
        
        remainingBoundariesMeshingEdge
        remainingBoundariesMeshingProperties
        
        % electrodeMeshingPropertiesAtRemainingDomain
        % remainingGeometryMeshingProperties
        
        meshingDomain
        domainMeshingProperties
        meshDistributionDDL % 1d
        meshDistributionRemaining % 1d
        
        % manualMesh
        testingMesh
        testingMeshBpeSurfaceMeshingEdge
        testingMeshBpeSurfaceMeshingProperties
        testingMeshTriangular
        testingMeshTriangularProperties
        
        %probes
        surfaceProbeTable
%         phiSurfaceProbe
%         cSurfaceProbe
        NSurfaceProbe
        iSurfaceProbe
        iCathodicSurfaceProbe
        iAnodicSurfaceProbe
        iTotalSurfaceProbe

        % 1d plots
        plotExporter1d
        
        standardPhiPlotGroup
        standardPhiPlot
        phiPBPlot
        standardConcentrationsPlotGroup
        standardConcentrationsPlot
        
        logConcentrationsPlotGroup
        logConcentrationsPlot
        phiAtSurfacePlotGroup
        phiAtSurfacePlot
        phiAtSurfacePBPlot

        logConcentrationsAtSurfacePlotGroup
        logConcentrationsAtSurfacePlot
        logPBConcentrationsAtSurfacePlot
        
        globalPlotGroup
        globalPlot
        
        multiPurpose1dPlotGroup
        multiPurpose1dPlot
        
        % 2d plots
        plotExporter2d  
        phiSurfacePlotGroup 
        cSurfacePlotGroup   
        phiContourPlotGroup 
        cContourPlotGroup   
        phiStreamlinePlotGroup 
        cStreamlinePlotGroup   
        phiArrowSurfacePlotGroup
        cArrowSurfacePlotGroup  
        phiArrowLinePlotGroup   
        cArrowLinePlotGroup    
        
        nStreamlinePlotGroup   
        nArrowSurfacePlotGroup  
        nArrowLinePlotGroup 
        
        iStreamlinePlotGroup   
        iArrowSurfacePlotGroup  
        iArrowLinePlotGroup 
        
        kappaSurfacePlotGroup
        kappaContourPlotGroup
        
        meshPlotGroup           
        
        % views
        standardView
        ddlView
        
        % datasets
        weResults
        ceResults
        bpeSurfaceResults
        zetaPlaneResults
        bulkBoundaryResults
        entireSurfaceResults
        centralCrossectionResults
        centralDDLCrossectionResults
        leftBpeEdgeCrossectionResults
        rightBpeEdgeCrossectionResults
        cathodeCrossectionResults
        anodeCrossectionResults
        halfDebyeLengthResults
        
        % studies and solvers
        stationaryStudy
        stationaryStudyStep1
        Stationary
        Parametric
        FullyCoupled
        
    end
    methods
        % get and set methods
        function v = get.RT(obj)
            v = jlh.Constants.MolarGasConstant*obj.T;
        end
        function v = get.UT(obj)
            v = obj.RT/obj.F;
        end
        function v = get.L(obj)
            v = obj.lambdaD/obj.epsilon;
        end
        function v = get.w(obj)
            v = obj.w_bulkLeft + obj.w_cathode + obj.w_insulatorLeft + ...
                obj.w_bpe + obj.w_insulatorRight + obj.w_anode + obj.w_bulkRight;
        end 
        function v = get.w_mesh(obj)
            v = obj.w_bpe/obj.nMeshChops;
        end
        function v = get.nMeshChops(obj)
            v = ceil(obj.w_bpe/obj.wMeshRatio/2)*2;
        end
        function v = get.nSegmentsPerChop(obj)
            v = round(obj.w_mesh/obj.elementSizeAtSurface);
        end
            
        function v = get.W(obj)
            v = obj.lambdaD*obj.w/obj.epsilon;
        end 
        function v = get.Wbpe(obj)
            v = obj.lambdaD*obj.w_bpe/obj.epsilon;
        end
        function v = get.WinsulatorLeft(obj)
            v = obj.lambdaD*obj.w_insulatorLeft/obj.epsilon;
        end
        function v = get.WinsulatorRight(obj)
            v = obj.lambdaD*obj.w_insulatorRight/obj.epsilon;
        end
        function v = get.Wcathode(obj)
            v = obj.lambdaD*obj.w_cathode/obj.epsilon;
        end
        function v = get.Wanode(obj)
            v = obj.lambdaD*obj.w_anode/obj.epsilon;
        end
        function v = get.WbulkLeft(obj)
            v = obj.lambdaD*obj.w_bulkLeft/obj.epsilon;
        end
        function v = get.WbulkRight(obj)
            v = obj.lambdaD*obj.w_bulkRight/obj.epsilon;
        end
%         function set.W(obj,v)
%             obj.w = obj.epsilon*v/obj.lambdaD;
%         end 
        function set.Wbpe(obj,v)
            obj.w_bpe = obj.epsilon*v/obj.lambdaD;
        end
        function set.WinsulatorLeft(obj,v)
            obj.w_insulatorLeft = obj.epsilon*v/obj.lambdaD;
        end
        function set.WinsulatorRight(obj,v)
            obj.w_insulatorRight = obj.epsilon*v/obj.lambdaD;
        end
        function set.Wcathode(obj,v)
            obj.w_cathode = obj.epsilon*v/obj.lambdaD;
        end
        function set.Wanode(obj,v)
            obj.w_anode = obj.epsilon*v/obj.lambdaD;
        end
        function set.WbulkLeft(obj,v)
            obj.w_bulkLeft = obj.epsilon*v/obj.lambdaD;
        end
        function set.WbulkRight(obj,v)
            obj.w_bulkRight = obj.epsilon*v/obj.lambdaD;
        end
        
        function v = get.PHI0(obj)
            v = obj.PHI_bpe;
        end
        function v = get.phi0(obj)
            v = obj.PHI_bpe / obj.UT;
        end
        function v = get.phi_bpe(obj)
            v = obj.PHI_bpe / obj.UT;
        end
        function v = get.delta_phi(obj)
            v = obj.deltaPHI / obj.UT;
        end
        function v = get.PHI_we(obj)
%             v = obj.PHI_bpe + obj.phiSymmetryFactor * obj.deltaPHI;
            v = obj.phiSymmetryFactor * obj.deltaPHI;
        end
        function v = get.PHI_ce(obj)
%             v = obj.PHI_bpe - (1-obj.phiSymmetryFactor) * obj.deltaPHI;
            v = - (1-obj.phiSymmetryFactor) * obj.deltaPHI;
        end
        function v = get.phi_we(obj)
%             v = obj.phi_bpe + obj.phiSymmetryFactor * obj.delta_phi;
            v = obj.phiSymmetryFactor * obj.delta_phi;
        end
        function v = get.phi_ce(obj)
%             v = obj.phi_bpe - (1-obj.phiSymmetryFactor) * obj.delta_phi;
            v = - (1-obj.phiSymmetryFactor) * obj.delta_phi;
        end
        function v = get.ionicStrength(obj)
            v = 1/2* sum(obj.z.^2.*obj.c_bulk); 
        end
        function v = get.lambdaD(obj)
            v = sqrt(jlh.Constants.VacuumPermittivity*obj.epsilon_r*obj.RT/(2*obj.F^2*obj.ionicStrength));
        end
        
        function v = get.c_ref(obj)
            v = obj.ionicStrength;
        end
        function v = get.D_ref(obj)
            v = sum(obj.D)/obj.numberOfSpecies;
        end
        function v = get.kappa(obj)
            v = obf.F^2/obj.RT*sum(obj.z.^2.*obj.c_bulk.^2.*obj.D);
        end
        
        function e = get.domainExpressions1d(obj)
            import jlh.*;
            import jlh.hf.*;
            logc = prepTerm('log(c_id)','c_id',obj.c_id);
            cx = prepTerm('c_idx','c_id',obj.c_id);
            C = prepTerm('c_id*c_ref','c_id',obj.c_id);
            logC = prepTerm('log(c_id*c_ref)','c_id',obj.c_id);
            Cx = prepTerm('(c_idx*c_ref/L)','c_id',obj.c_id);
            e = {     {'phi'}, {'phix'},...
                        {'phi*UT'}, {'phix*UT/L'},...
                        obj.c_id, logc, cx, ...
                        C,logC,Cx,...
                        obj.nx_id,{'ix'},{'kappa_loc'} };
        end
        function e = get.surfaceExpressions1d(obj)
            e = { {'i_total'},...
                {'i_cathodic', 'i_anodic'},...
                {'log(abs(i_cathodic))','log(abs(i_anodic))'},...
                obj.i_id, obj.i_dimless_id, ...
                obj.N_id, obj.N_dimless_id};
        end
        
        function t = e2t(~,e) % expression to text
            import jlh.hf.*;
            t = flattenCell(e);
            t = strrep(t,'(','_');
            t = strrep(t,')','');
            t = strrep(t,'/','_by_');
            t = strrep(t,'*','_times_');
            t = strrep(t,'-','_minus_');
            t = strrep(t,'+','_plus_');
        end
        
        function updateFluxEvaluations1d(obj,dset)
            obj.updateEvaluations1d(dset,[obj.N_dimless_id obj.nx_id]);
        end
        
        function n = get.nMultiPurpose1dPlots(obj)
            n  = max(obj.numberOfSpecies,obj.nReactions);
        end

        function dset = getLatestDataset(obj)
            dset = strtrim(jlh.hf.lastRow(char(obj.m.result.dataset.tags)));
        end
        function f = makeFileName(obj,n)
            f = [obj.projectPath,'/',n];
        end
        
        function ap = getIdentityPairsForComponent(obj, component) 
            pairTags = obj.m.pair.tags;
            ap = {};
            for j=1:numel(pairTags)
                c = char(obj.m.pair(pairTags(j)).model);
                if strcmp(c,component)
                    ap = [ap, char(pairTags(j))];
                end
            end
        end
    
        obj = prepareIdentifiers(obj)
        obj = createModel(obj)
        obj = createGeometry(obj)
        obj = createChoppedGeometry(obj)
        obj = createSelections(obj)
        obj = createPhysics(obj)
        
        create0dComponent(obj)
        
        obj = makeChoppedGeometry(obj)
        obj = makeSimpleGeometry(obj)
        
        obj = updateModel(obj)
        
        obj = updateParameters(obj)
        obj = updateFunctions(obj)
        obj = updateGeometry(obj)
        obj = updateChoppedGeometry(obj)
        obj = updateSelections(obj)
        obj = updateOperators(obj)
        obj = updatePhysics(obj)
        obj = updateMesh(obj)
        obj = updateMeshMapped(obj)

        obj = createMappedMesh(obj);
        obj = updateMappedMesh(obj);
        obj = finalizeMappedMesh(obj);
        
        obj = createMappedMeshForSimpleGeometry(obj);
        obj = updateMappedMeshForSimpleGeometry(obj);

        
        obj = updateChoppedMesh(obj)
        obj = replicateMeshPrototype(obj)
        obj = finalizeChoppedMesh(obj)
        obj = updateOperatorsForChoppedGeometry(obj)
        obj = updatePhysicsForChoppedGeometry(obj)
        obj = updateDatasetsForChoppedGeometry(obj)
        
        obj = createPhysicsForWeakForm(obj)
        obj = updatePhysicsForWeakForm(obj)
    
        obj = updateViews(obj)
        obj = updateDatasets(obj,dset)
        
        obj = addParametricStudy(obj,parameterName,parameterValues)
        obj = addBpeStudy(obj,parameterName,parameterValues)
        obj = addBpeStudyStepSequence(obj,parameterName,parameterValues,disablePhysics)

        obj = mesh1D(obj)

        create1dFluxModel(obj)
        create1dFluxGeometry(obj)
        create1dFluxSelections(obj)
        create1dFluxOperators(obj)
        create1dFluxPhysics(obj)
        
        update1dFluxGeometry(obj)
        update1dFluxSelections(obj)
        update1dFluxOperators(obj)
        update1dFluxPhysics(obj)
        
        obj = sweepChargeDensityCoupling(obj)

        obj = newProject(obj,name)
        obj = clearProject(obj)
        obj = deleteProject(obj)
        obj = loadMph(obj)
        obj = loadState(obj)
        obj = load1dState(obj)
        obj = saveState(obj)
        
        obj = sampleProject(obj)
        obj = initializeSampleProject(obj)
        obj = runSampleProject(obj)
        
        obj = setKrawiecCaseGeneralParameters(obj)
        obj = setKrawiecCaseReactionParameters(obj)
        obj = setZhangCaseParameters(obj)
        obj = setDuval2001BipolarCaseParameters(obj,c)
        obj = setDhopeshwarkar2008ElectrokineticsCaseParameters(obj)
        obj = setSampleCaseParameters(obj)

        obj = setFeederElectrodeCurrentFluxAndOpenBoundaryConditions(obj)
        obj = setFeederElectrodeCurrentFluxAndBPEPotentialBC(obj)
        obj = setSurfaceOverpotentialExpression(obj)
        
        activateBpe(obj)
        disableBpe(obj)
        
        % plotting
        obj = initializeExporter(obj)
        obj = update1dPlots(obj)
        obj = update1dParametricPlots(obj,dset)
        obj = update2dPlots(obj,dset)
        obj = update2dParametricPlots(obj,dset)
        obj = updateGlobalPlots(obj,dset,par)
        obj = updateGlobalPlotsBySolnum(obj,dset,par,plotDset);
        obj = updatePlotsDimensionless(obj,dset)
        obj = plotCurrents(obj)
        obj = plotAlongCrossection(obj,dset,dsetFolder)
        obj = plotPresetCrossections(obj,dset)
        obj = sweepHorizontalCrossection(obj,dset,stepSize,left,right,lower,upper)
%         obj = sweepHorizontalCrossectionAtBpe(obj,dset,stepSize)
        obj = sweepVerticalCrossection(obj,dset,stepSize,left,right,lower,upper)
        obj = iterateStandardPlots(obj)
        
        obj = createEvaluations1d(obj);
        obj = updateEvaluations1d(obj,dset,eval,table_id);
        obj = plot1dSolution(obj,dset,dsetFolder);
        obj = plot1dSolutionSurfaceParametric(obj,dset,par,dsetFolder);
        obj = evaluate1dSolution(obj,dset,dsetFolder);
        txtFileName = exportSolution1d(obj,dset,dsetFolder);
        
        function obj = savePlot(obj,f,fileName,nRows,nCols)
            if isempty(obj.projectPath)
                fileName = ['img/',fileName];
            else
                fileName = obj.makeFileName(fileName);
            end
            set(f,'PaperPosition',5*[0 0 nCols nRows/obj.widthToHeight]);
            print(f,'-r300','-dpng',fileName);
        end     
    end
end