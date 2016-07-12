% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

caseStudyParameterFile = 'parameters_duval2003coupling.m';
caseStudyTitleSuffix = '_duval2003coupling';

makeGeometryPartFileMacro

runTertiaryCurrentDistribution2dWithNetCurrentConstraint
runDilutedSpeciesAndElectrostatics1dSweep
runDilutedSpeciesAndElectrostatics2dPureDirichletBC
runDilutedSpeciesAndElectrostatics2d