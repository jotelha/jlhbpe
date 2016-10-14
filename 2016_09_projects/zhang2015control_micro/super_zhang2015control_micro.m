caseStudyParameterFile = 'parameters_zhang2015control_micro.m';

caseStudyTitleSuffix = '_zhang2015control_micro';

hMaxFactor = 0.001;
makeGeometryPartFile

hMaxFactor = 0.01;
makeGeometryPartFile

% meshFile = 'copySimpleBulkGeometryMeshCoarse_case_c.m';
% meshFile = 'copySimpleBulkGeometryRefinedMesh.m';

pname = { 'sweep' 'DeltaPHI' };
plistarr = { 'range(1,1,3)' '0.1563,0.3166,0.5641' };

batchTertiaryCurrentDistribution2d

batchTertiaryCurrentDistribution2dWithNetCurrentConstraint

hMaxFactor = 0.001;
batchPDE1dSweep

batchDilutetSpeciesAndElectrostatics2dPureDirichlet

batchDilutetSpeciesAndElectrostatics2d
% runDilutedSpeciesAndElectrostatics1dSweep
% runDilutedSpeciesAndElectrostatics2dPureDirichletBC
% runDilutedSpeciesAndElectrostatics2d