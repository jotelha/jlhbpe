% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

% meshing by hand

caseStudyParameterFile = 'parameters_duval2003faradaic_micro2.m';

c = 4;
caseStudyTitleSuffix = '_duval2003faradaic_micro2_case_d';

hMaxFactor = 0.1;
makeGeometryPartFile

% meshFile = 'copySimpleBulkGeometryMeshCoarse_case_c.m';
% meshFile = 'copySimpleBulkGeometryRefinedMesh.m';

pname = { 'sweep' 'DeltaPHI' };
plistarr = { 'range(1,1,1)' 'range(0.2,0.2,0.2)' };

batchTertiaryCurrentDistribution2d

batchTertiaryCurrentDistribution2dWithNetCurrentConstraint

hMaxFactor = 0.01;
batchPDE1dSweep
% runDilutedSpeciesAndElectrostatics1dSweep

%batchDilutetSpeciesAndElectrostatics2dPureDirichlet

batchDilutetSpeciesAndElectrostatics2d