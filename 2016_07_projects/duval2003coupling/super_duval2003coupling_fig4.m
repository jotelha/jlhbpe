% caseStudyParameterFile = 'parameters_duval2001bipolar.m';
% caseStudyTitleSuffix = '_duval2001bipolar';

caseStudyParameterFile = 'parameters_duval2003coupling_fig4.m';
caseStudyTitleSuffix = '_duval2003coupling';

makeGeometryPartFileMacro

% runTertiaryCurrentDistribution2d
% runTertiaryCurrentDistribution2dWithNetCurrentConstraint
% runDilutedSpeciesAndElectrostatics1dSweep
% runDilutedSpeciesAndElectrostatics2dPureDirichletBC
% runDilutedSpeciesAndElectrostatics2d

pname = { 'sweep' 'DeltaPHI' };
plistarr = { 'range(1,1,10)' 'range(0.5,0.5,5)' };

batchTertiaryCurrentDistribution2d



batchTertiaryCurrentDistribution2dWithNetCurrentConstraint