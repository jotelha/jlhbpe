% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

caseStudyParameterFile = 'parameters_duval2003faradaic.m';

c = 1;
caseStudyTitleSuffix = '_duval2003faradaic_case_a';


makeGeometryPartFileMacro

batchTertiaryCurrentDistribution2d
% runDilutedSpeciesAndElectrostatics1dSweep
% runDilutedSpeciesAndElectrostatics2dPureDirichletBC
% runDilutedSpeciesAndElectrostatics2d