function obj = load1dState(obj)
    %% prepare paths
    % projectName and projectFolderName have to be set

    if isempty(obj.projectName) %default: latest project, tag in workspace
        projects = dir('dat\*');
        obj.projectName = projects(end).name;
%     elseif exist('model_tag','var') ~= 1 % if projectName is given, but no model_tag
%         model_tag = projectName;
    end
    if isempty(obj.projectPath) %default: folder and project have same names
        obj.projectPath = ['dat\',obj.projectName];
    end
    if isempty(obj.model_tag) %default: latest mph file in project folder
        mphFiles= dir(sprintf('%s\\*.mph',obj.projectPath));
        mphFile = [obj.projectPath,'\',mphFiles(end).name];
        obj.model_tag = strrep(mphFiles(end).name,'.mph','');
    else
        mphFile = [obj.projectPath,'\',obj.model_tag,'.mph'];
    end
    %projectName =  [datestr(date,'yyyy_mm_dd_'),model_tag];

    fprintf('Loading state from folder %s:\n', obj.projectPath);
%     projectPath = ['dat/',projectFolderName];
%     mphFile = [projectPath,'/',projectName,'.mph'];

    matFile = [obj.projectPath,'\',obj.projectName,'.mat']; % has matlab variables
%     datFile = [projectPath,'/',projectName,'_dat.mat']; % has solution data

%     generalParameterFile = [projectPath,'/setGeneralParameters.m'];
%     reactionParameterFile = [projectPath,'/setReactionParameters.m'];

    %% set parameters
%     fprintf('Loading general parameters from %s...\n',generalParameterFile);
%     run(generalParameterFile);
% 
%     fprintf('Loading reaction parameters from %s...\n',reactionParameterFile);
%     run(reactionParameterFile);
% 
%     fprintf('Rebuilding identifieres for COMSOL-internal use...\n');
%     prepareIdentifiers
% 
    fprintf('Loading MATLAB class from %s...\n',matFile);
    dat = load(matFile);

%     %% load solution data
%     fprintf('Loading solution data extracted from COMSOL earlier from %s...\n',matFile);
%     load(datFile);

%     %% load COMSOL model
%     if exist('model_tag','var')
%         fprintf('Loading COMSOL model %s to server, assigning tag "%s"...\n',mphFile,model_tag);
%         m = mphload(mphFile,model_tag);
%     else
%         fprintf('Loading COMSOL model %s to server assigning default tag"...\n',mphFile);
%         m = mphload(mphFile);
%     end
    fprintf('Loading COMSOL model %s to server, assigning tag "%s"...\n',mphFile,obj.model_tag);
    m = mphload(mphFile,obj.model_tag);
    obj = dat.obj;
    obj.m = m;
        
    %% restore functions
    obj.phiInterpolation = m.func('phiInterpolation'); 
    
    obj.phi_pb = m.func('phi_pb');
    obj.phi_pbx = m.func('phi_pbx');
    obj.c_pb = cell(obj.numberOfSpecies,1);
    obj.cInterpolation = cell(obj.numberOfSpecies,1);
    for i = 1:obj.numberOfSpecies   
        obj.c_pb{i} = m.func(obj.c_pb_id{i});
        obj.cInterpolation{i} = m.func(obj.cInterpolation_id{i});
    end

    %% restore geometry
    obj.geom = m.geom('geom');
%     space = geom.feature('space');
    obj.ddl             = obj.geom.feature('ddl');
    obj.space           = obj.geom.feature('space');
    obj.leftEdgeOfBpe   = obj.geom.feature('leftEdgeOfBpe');
    obj.rightEdgeOfBpe  = obj.geom.feature('rightEdgeOfBpe');
%     obj.leftEdgeOfCathode  = obj.geom.feature('leftEdgeOfCathode');
%     obj.rightEdgeOfCathode  = obj.geom.feature('rightEdgeOfCathode');
%     obj.leftEdgeOfAnode  = obj.geom.feature('leftEdgeOfAnode');
%     obj.rightEdgeOfAnode  = obj.geom.feature('rightEdgeOfAnode');

    %% restore selections
%     surfaceNode = m.selection('surfaceNode');
%     zetaNode = m.selection('zetaNode');
%     bulkNode = m.selection('bulkNode');

%     obj.leftBoundaryOfSurface = obj.m.selection('leftBoundaryOfSurface');
%     obj.rightBoundaryOfSurface = obj.m.selection('rightBoundaryOfSurface');
%     obj.leftBoundaryOfZetaPlane = obj.m.selection('leftBoundaryOfZetaPlane');
%     obj.rightBoundaryOfZetaPlane = obj.m.selection('rightBoundaryOfZetaPlane');
    
    % edge selections
    obj.allBoundaries                   = obj.m.selection('allBoundaries');
    obj.bpeSurface                      = obj.m.selection('bpeSurface');
    obj.bulkBoundary                    = obj.m.selection('bulkBoundary');
    obj.zetaPlane                       = obj.m.selection('zetaPlane'); % edge at some distance from BPE (one Debye length), representing the zeta plane
%     obj.lateralBoundaryAtFirstDebyeLength    = obj.m.selection('lateralBoundaryAtFirstDebyeLength'); 
%     obj.lateralBoundaryAtRemainingDomain     = obj.m.selection('lateralBoundaryAtRemainingDomain');
%     obj.lateralBoundary                 = obj.m.selection('lateralBoundary');
%     obj.upperBoundary                   = obj.m.selection('upperBoundary');
%     obj.electrodes                      = obj.m.selection('electrodes');
%     obj.workingElectrode                = obj.m.selection('workingElectrode');
%     obj.counterElectrode                = obj.m.selection('counterElectrode');
%     obj.insulator                       = obj.m.selection('insulator');
%     obj.insulatorAdjacentToBpe          = obj.m.selection('insulatorAdjacentToBpe');
    obj.entireSurface                   = obj.m.selection('entireSurface');
    obj.reactingSurface                 = obj.m.selection('reactingSurface');
    obj.regionOfFirstDebyeLength        = obj.m.selection('regionOfFirstDebyeLength');
    obj.regionRemaining                 = obj.m.selection('regionRemaining');

    %% restore local definitions
%     projectZetaPotential = m.cpl('projectZetaPotential');
%     projectSurfaceToBulkExit = m.cpl('projectSurfaceToBulkExit');
    obj.projectReactionPlaneToSurface = obj.m.cpl('projectReactionPlaneToSurface');

    obj.integrateSurface    = obj.m.cpl('integrateSurface');
    obj.integrateBulkExit   = obj.m.cpl('integrateBulkExit');
    obj.maximumOnDomain     = obj.m.cpl('maximumOnDomain');

    % variables
    obj.domainVariables     = obj.m.variable('domainVariables');
    obj.surfaceVariables    = obj.m.variable('surfaceVariables');
    obj.bulkVariables    = obj.m.variable('bulkVariables');
    obj.electrodeVariables  = obj.m.variable('electrodeVariables');
    obj.weVariables         = obj.m.variable('weVariables');
    obj.ceVariables         = obj.m.variable('ceVariables');

    %% restore physics
    obj.NernstPlanckEquation    = obj.m.physics('NernstPlanckEquation');
    obj.FluxAtSurface           = obj.NernstPlanckEquation.feature('FluxAtSurface');
    obj.FluxAtElectrodes   = obj.NernstPlanckEquation.feature('FluxAtElectrodes');
    obj.FluxAtBulkBoundary          = obj.NernstPlanckEquation.feature('FluxAtBulkBoundary');
    obj.BulkConcentrationsBC    = obj.NernstPlanckEquation.feature('BulkConcentrationsBC');
    
    obj.PoissonEquation         = obj.m.physics('PoissonEquation');
    obj.bpeChargeDensityBC      = obj.PoissonEquation.feature('bpeChargeDensityBC');
    obj.bpePotentialBC          = obj.PoissonEquation.feature('bpePotentialBC');
    obj.bulkPotentialBC         = obj.PoissonEquation.feature('bulkPotentialBC');
    obj.electrodePotentialBC    = obj.PoissonEquation.feature('electrodePotentialBC');  
    obj.bulkBoundaryPotentialGradientBC    = obj.PoissonEquation.feature('bulkBoundaryPotentialGradientBC');  
    obj.entireSurfacePotentialGradientBC   = obj.PoissonEquation.feature('entireSurfacePotentialGradientBC');  
    obj.insulatorSurfacePotentialGradientBC   = obj.PoissonEquation.feature('insulatorSurfacePotentialGradientBC');  
    
    obj.zeroSurfaceCurrent          = obj.m.physics('zeroSurfaceCurrent');
    obj.zeroSurfaceCurrentEquation  = obj.zeroSurfaceCurrent.feature('zeroSurfaceCurrentEquation');

    %% restore meshes
    obj.standardMesh                        = obj.m.mesh('standardMesh');    
    obj.meshingEdge                         = obj.standardMesh.feature('meshingEdge'); % for 1d model
    obj.meshDistributionDDL                 = obj.meshingEdge.feature('meshDistributionDDL');
    obj.meshDistributionRemaining           = obj.meshingEdge.feature('meshDistributionRemaining');
%     obj.bpeSurfaceMeshingEdge               = obj.standardMesh.feature('bpeSurfaceMeshingEdge');
%     obj.bpeSurfaceMeshingProperties         = obj.bpeSurfaceMeshingEdge.feature('bpeSurfaceMeshingProperties');
% 
%     obj.zetaPlaneMeshingEdge                = obj.standardMesh.feature('zetaPlaneMeshingEdge');
%     obj.zetaPlaneMeshingProperties          = obj.zetaPlaneMeshingEdge.feature('zetaPlaneMeshingProperties');
%     
%     obj.electrodeMeshingEdgeAtFirstDebyeLength          = obj.standardMesh.feature('electrodeMeshingEdgeAtFirstDebyeLength');
%     obj.electrodeMeshingPropertiesAtFirstDebyeLength    = obj.electrodeMeshingEdgeAtFirstDebyeLength.feature('electrodeMeshingPropertiesAtFirstDebyeLength');

%     obj.remainingBoundariesMeshingEdge               = obj.standardMesh.feature('remainingBoundariesMeshingEdge');
%     obj.remainingBoundariesMeshingProperties         = obj.remainingBoundariesMeshingEdge.feature('remainingBoundariesMeshingProperties');

%     obj.meshingDomain                       = obj.standardMesh.feature('meshingDomain'); 
%     obj.domainMeshingProperties             = obj.meshingDomain.feature('domainMeshingProperties');
    
    % a testing mesh
%     obj.testingMesh                             = obj.m.mesh('manualMesh');
%     obj.testingMeshBpeSurfaceMeshingEdge        = obj.testingMesh.feature('testingMeshBpeSurfaceMeshingEdge');
%     obj.testingMeshBpeSurfaceMeshingProperties  = obj.testingMeshBpeSurfaceMeshingEdge.feature('testingMeshBpeSurfaceMeshingProperties');
% 
%     obj.testingMeshTriangular               = obj.testingMesh.feature('testingMeshTriangular'); 
%     obj.testingMeshTriangularProperties     = obj.testingMeshTriangular.feature('testingMeshTriangularProperties');
% 

    %% restore last study and solution
%     studyTags = m.study.tags;
%     if numel(studyTags) > 0
%         stationaryStudy = m.study(studyTags(end));
%         stationaryStudyStep1 = stationaryStudy.feature('stationaryStudyStep1');
%     end
% 
%     solTags = m.sol.tags;
%     if numel(solTags) > 0
%         sol = m.sol(solTags(end));
%         Stationary = sol.feature('s1');
%         FullyCoupled = Stationary.feature('fc1');
%     end
 %% restore probes
%     obj.surfaceProbeTable = obj.m.result.table('probeTable');
% 
%     obj.phiSurfaceProbe = obj.m.result.numerical('phiSurfaceProbe');
% 
%     obj.cSurfaceProbe = cell(obj.numberOfSpecies,1);
%     obj.NSurfaceProbe = cell(obj.numberOfSpecies,1);
%     for i=1:obj.numberOfSpecies  
%         obj.cSurfaceProbe{i} = obj.m.result.numerical(obj.cSurfaceProbe_id{i});
%         obj.NSurfaceProbe{i} = obj.m.result.numerical(obj.NSurfaceProbe_id{i});
%     end
% 
%     obj.iCathodicSurfaceProbe   = obj.m.result.numerical('iCathodic');
%     obj.iAnodicSurfaceProbe     = obj.m.result.numerical('iAnodic');
%     obj.iTotalSurfaceProbe      = obj.m.result.numerical('iTotal');

    %% restore 1d plots
    
    obj.plotExporter1d          = obj.m.result.export('plotExporter1d');

    obj.standardPhiPlotGroup    = obj.m.result('standardPhiPlotGroup');
    obj.standardPhiPlot         = obj.standardPhiPlotGroup.feature('standardPhiPlot');
    obj.phiPBPlot               = obj.standardPhiPlotGroup.feature('phiPBPlot');


    obj.standardConcentrationsPlotGroup     = obj.m.result('standardConcentrationsPlotGroup');
    obj.standardConcentrationsPlot          = cell(1,obj.numberOfSpecies);
    for i=1:obj.numberOfSpecies
        obj.standardConcentrationsPlot{i} = obj.standardConcentrationsPlotGroup.feature(obj.standardConcentrationsPlot_id{i});
    end

    obj.logConcentrationsPlotGroup = obj.m.result('logConcentrationsPlotGroup');
    obj.logConcentrationsPlot = cell(1,obj.numberOfSpecies);
    for i=1:obj.numberOfSpecies
        obj.logConcentrationsPlot{i} = obj.logConcentrationsPlotGroup.feature(obj.logConcentrationsPlot_id{i});
    end

    obj.globalPlotGroup    = obj.m.result('globalPlotGroup');
    obj.globalPlot         = obj.globalPlotGroup.feature('globalPlot');
%     %% restore obj.ddlView. plots of ddl region just at surface, phi
%     % potential
%     obj.phiAtSurfacePlotGroup = obj.m.result.create('phiAtSurfacePlotGroup',1);
%     obj.phiAtSurfacePlot = obj.phiAtSurfacePlotGroup.create('phiAtSurfacePlot','LineGraph');
%     obj.phiAtSurfacePBPlot = obj.phiAtSurfacePlotGroup.feature.duplicate('phiAtSurfacePBPlot','phiAtSurfacePlot');
% 
%     %% create plots of ddl region just at surface, log concentrations
%     obj.logConcentrationsAtSurfacePlotGroup = obj.m.result.create('logConcentrationsAtSurfacePlotGroup',1);
%     obj.logConcentrationsAtSurfacePlot = cell(1,obj.numberOfSpecies);
%     obj.logConcentrationsAtSurfacePlot_id = nprintf('logConcentrationsAtSurfacePlot%d',obj.numberOfSpecies);
%     obj.logPBConcentrationsAtSurfacePlot = cell(1,obj.numberOfSpecies);
%     obj.logPBConcentrationsAtSurfacePlot_id = nprintf('logPBConcentrationsAtSurfacePlot%d',obj.numberOfSpecies);
%     for i=1:obj.numberOfSpecies
%         obj.logConcentrationsAtSurfacePlot{i} = obj.logConcentrationsAtSurfacePlotGroup.create(obj.logConcentrationsAtSurfacePlot_id{i},'LineGraph');
%         obj.logPBConcentrationsAtSurfacePlot{i} = obj.logConcentrationsAtSurfacePlotGroup.create(obj.logPBConcentrationsAtSurfacePlot_id{i},'LineGraph');
%     end

    %% restored 2d plots
%     obj.plotExporter2d      = obj.m.result.export('plotExporter2d');
%     
%     obj.phiSurfacePlotGroup = obj.m.result('phiSurfacePlotGroup');
%     obj.cSurfacePlotGroup   = obj.m.result('cSurfacePlotGroup');
%     
%     obj.phiContourPlotGroup = obj.m.result('phiContourPlotGroup');
%     obj.cContourPlotGroup   = obj.m.result('cContourPlotGroup');
%    
%     obj.phiStreamlinePlotGroup  = obj.m.result('phiStreamlinePlotGroup');
%     obj.cStreamlinePlotGroup    = obj.m.result('cStreamlinePlotGroup');
% 
%     obj.phiArrowSurfacePlotGroup= obj.m.result('phiArrowSurfacePlotGroup');
%     obj.cArrowSurfacePlotGroup  = obj.m.result('cArrowSurfacePlotGroup');
%     
%     obj.phiArrowLinePlotGroup   = obj.m.result('phiArrowLinePlotGroup');
%     obj.cArrowLinePlotGroup     = obj.m.result('cArrowLinePlotGroup');
%     
%     obj.nStreamlinePlotGroup    = obj.m.result('nStreamlinePlotGroup');
%     obj.nArrowSurfacePlotGroup  = obj.m.result('nArrowSurfacePlotGroup');
%     obj.nArrowLinePlotGroup     = obj.m.result('nArrowLinePlotGroup');
%  
%     obj.iStreamlinePlotGroup    = obj.m.result('iStreamlinePlotGroup');
%     obj.iArrowSurfacePlotGroup  = obj.m.result('iArrowSurfacePlotGroup');
%     obj.iArrowLinePlotGroup     = obj.m.result('iArrowLinePlotGroup');
%     
%     obj.kappaSurfacePlotGroup = obj.m.result('kappaSurfacePlotGroup');    
%     obj.kappaContourPlotGroup = obj.m.result('kappaContourPlotGroup');
% 
%     obj.meshPlotGroup           = obj.m.result('meshPlotGroup');
    
    %% restore views
%     obj.standardView = obj.m.view('standardView');
%     obj.ddlView      = obj.m.view('ddlView');
   
    %% restore datasets
%     obj.weResults                   = obj.m.result.dataset('weResults');
%     obj.ceResults                   = obj.m.result.dataset('ceResults');
%     obj.bpeSurfaceResults           = obj.m.result.dataset('bpeSurfaceResults');
%     obj.zetaPlaneResults            = obj.m.result.dataset('zetaPlaneResults');
%     obj.bulkBoundaryResults         = obj.m.result.dataset('bulkBoundaryResults');    
%     obj.entireSurfaceResults        = obj.m.result.dataset('entireSurfaceResults');
% 	obj.centralCrossectionResults   = obj.m.result.dataset('centralCrossectionResults');
% 	obj.centralDDLCrossectionResults= obj.m.result.dataset('centralDDLCrossectionResults');
% 	
%     obj.leftBpeEdgeCrossectionResults   = obj.m.result.dataset('leftBpeEdgeCrossectionResults');
% 	obj.rightBpeEdgeCrossectionResults  = obj.m.result.dataset('rightBpeEdgeCrossectionResults');
% 	obj.cathodeCrossectionResults       = obj.m.result.dataset('cathodeCrossectionResults');
% 	obj.anodeCrossectionResults         = obj.m.result.dataset('anodeCrossectionResults');
% 	obj.halfDebyeLengthResults          = obj.m.result.dataset('halfDebyeLengthResults');
end