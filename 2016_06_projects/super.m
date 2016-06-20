% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

caseStudyParameterFile = 'parameters_duval2003faradaic.m';
caseStudyTitleSuffix = '_duval2003faradaic';

makeGeometryPartFile 

runTertiaryCurrentDistribution2d
runDilutedSpeciesAndElectrostatics1dSweep
runDilutedSpeciesAndElectrostatics2d