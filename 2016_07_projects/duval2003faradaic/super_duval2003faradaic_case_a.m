% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

caseStudyParameterFile = 'parameters_duval2003faradaic.m';

c = 1;
caseStudyTitleSuffix = '_duval2003faradaic_case_a';


makeGeometryPartFileMacro

meshFile = 'copySimpleBulkGeometryMeshCoarse_case_a.m';

pname = { 'sweep' 'DeltaPHI' };
plistarr = { 'range(1,1,10)' 'range(0.2,0.2,2)' };

batchTertiaryCurrentDistribution2d

batchTertiaryCurrentDistribution2dWithNetCurrentConstraint

hMaxFactor = 0.001;
batchPDE1dSweep
% runDilutedSpeciesAndElectrostatics1dSweep
% runDilutedSpeciesAndElectrostatics2dPureDirichletBC
% runDilutedSpeciesAndElectrostatics2d