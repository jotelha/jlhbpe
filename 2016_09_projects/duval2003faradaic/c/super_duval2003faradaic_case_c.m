% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

% meshing by hand

caseStudyParameterFile = 'parameters_duval2003faradaic.m';

c = 3;
caseStudyTitleSuffix = '_duval2003faradaic_case_c';


makeGeometryPartFileMacro

% meshFile = 'copySimpleBulkGeometryMeshCoarse_case_c.m';
% meshFile = 'copySimpleBulkGeometryRefinedMesh.m';

pname = { 'sweep' 'DeltaPHI' };
plistarr = { 'range(1,1,10)' 'range(0.2,0.2,2)' };

batchTertiaryCurrentDistribution2d

batchTertiaryCurrentDistribution2dWithNetCurrentConstraint

hMaxFactor = 0.001;
batchPDE1dSweep

% runDilutedSpeciesAndElectrostatics1dSweep
% runDilutedSpeciesAndElectrostatics2dPureDirichletBC
% runDilutedSpeciesAndElectrostatics2d