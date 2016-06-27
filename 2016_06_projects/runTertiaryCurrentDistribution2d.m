import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end

caseStudyTitle = 'tertiaryCurrentdDistribution2d';
createEmptyProject;

% or, to load an existing project from server or from file
% tag = 'xyz'
% loadExistingProject

% prepare text files
% makeParameterFile
% makeVariablesFiles

%% load parameters, geometry and mesh parts
% T = readtable('geometryPartsFile.txt','ReadRowNames',true,'ReadVariableNames',false,'Delimiter',' ');
T = readtable('globalFiles.txt','ReadRowNames',false,'ReadVariableNames',false,'Delimiter',' ');
% writetable(T,'geometryPartsFile.txt','WriteRowNames',true,'WriteVariableNames',false,'Delimiter',' ');
C = table2cell(T);
files = containers.Map(C(:,1),C(:,2));
% geometryPartsMphFile = T.(1)('geometryPartsMphFile');
% simpleBulkGeometryRefinedMeshFile = T.(1)('simpleBulkGeometryRefinedMeshFile');
% simpleDdlGeometryRefinedMeshFile = T.(1)('simpleDdlGeometryRefinedMeshFile');

model = m.m;

model.param().loadFile(files('parameterFile'));

% mesh parts

model.modelNode.create('mcomp1', 'MeshComponent');
model.geom.create('mgeom1', 2); % mesh geometry
model.mesh.create('simpleBulkGeometryRefinedMeshPart', 'mgeom1');
model.mesh('simpleBulkGeometryRefinedMeshPart').create('imp1', 'Import');
model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('filename', files('simpleBulkGeometryRefinedMeshFile'));
model.mesh('simpleBulkGeometryRefinedMeshPart').run;
% model.mesh('mpart1').feature('imp1').set('facepartition', 'auto');

% model.modelNode.create('mcomp2', 'MeshComponent');
% model.geom.create('mgeom2', 2); % mesh geometry
% model.mesh.create('simpleDdlGeometryRefinedMeshPart', 'mgeom2');
% model.mesh('simpleDdlGeometryRefinedMeshPart').create('imp1', 'Import');
% model.mesh('simpleDdlGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
% model.mesh('simpleDdlGeometryRefinedMeshPart').feature('imp1').set('filename', files('simpleDdlGeometryRefinedMeshFile'));
% model.mesh('simpleDdlGeometryRefinedMeshPart').run;

%%
makeTertiaryCurrentDistribution2dComponent

makeSimpleStudy
% tags of created features in sol_id, study_id, studyStep_id, compile_id,
% variables_id, solver_id and store_id
model.sol('sol1').feature('s1').feature('dDef').set('ooc', 'on');


model.sol(sol_id).run(store_id);

storedSolutions = model.study(study_id).getSolverSequences('Stored');
solInfo = mphsolinfo(model,'soltag',char(storedSolutions(1)));
dset = char(solInfo.dataset);

m.saveState;

createPlots

%%
% set up expressions to plot
% small letters: dimensionles, capital letters: dimensional

logc = prepTerm('log(c_id)','c_id',m.c_id);
% cx = prepTerm('c_idx','c_id',obj.c_id);
% cy = prepTerm('c_idy','c_id',obj.c_id);
cx = prepTerm('tcdee.grad_c_idx','c_id',m.c_id);
cy = prepTerm('tcdee.grad_c_idy','c_id',m.c_id);
absc  = prepTerm('sqrt(cx^2+cy^2)','cx','cy',cx,cy);
% dimensional
C = prepTerm('c_id/c_ref','c_id',m.c_id);
logC = prepTerm('log(c_id/c_ref)','c_id',m.c_id);
Cx = prepTerm('(tcdee.grad_c_idx/c_ref*L)','c_id',m.c_id);
Cy = prepTerm('(tcdee.grad_c_idy/c_ref*L)','c_id',m.c_id);
absC  = prepTerm('sqrt(Cx^2+Cy^2)','Cx','Cy',Cx,Cy);

absn = prepTerm('sqrt(nx_id^2+ny_id^2)','nx_id','ny_id',m.nx_id,m.ny_id);
absN = prepTerm('sqrt(Nx_id^2+Ny_id^2)','Nx_id','Ny_id',m.Nx_id,m.Ny_id);

absDiffusiveFlux = prepTerm('tcdee.dfluxMag_c_id','c_id',m.c_id);
diffusiveFluxX = prepTerm('tcdee.dflux_c_idx','c_id',m.c_id);
diffusiveFluxY = prepTerm('tcdee.dflux_c_idy','c_id',m.c_id);
electrophoreticFluxX = prepTerm('tcdee.mflux_c_idx','c_id',m.c_id);
electrophoreticFluxY = prepTerm('tcdee.mflux_c_idy','c_id',m.c_id);
absElectrophoreticFlux = prepTerm('sqrt(efx^2+efy^2)','efx','efy',electrophoreticFluxX,electrophoreticFluxY);
absTotalFlux = prepTerm('tcdee.tfluxMag_c_id','c_id',m.c_id);
totalFluxX = prepTerm('tcdee.tflux_c_idx','c_id',m.c_id);
totalFluxY = prepTerm('tcdee.tflux_c_idy','c_id',m.c_id);
       
            
plots = containers.Map;
% plots(title) = { {expression1, expression2, ...}, ylabel };
plots('phi')    = { {'phi/UT'}, 'phi / U_T'};
plots('phix')   = { {'phix/UT*L'}, 'phi_x * L / U_T'};
plots('phiy')   = { {'phiy/UT*L'}, 'phi_y * L / U_T'};
plots('_PHI')   = { {'phi'}, 'phi / V' };
plots('_PHIx')  = { {'phix'}, 'phi_x / V m^{-1}'};
plots('_PHIy')  = { {'phiy'}, 'phi_y / V m^{-1}'};

plots('c')      = { C, 'c / c_ref' };
plots('logc')   = { logC, 'log(c/c_ref)'};
plots('absc')   = { absC,'|grad(c)| * L / c_ref'}; 
plots('cx')     = { Cx,'c_x * L / c_ref'};
plots('cy')     = { Cy, 'c_y * L / c_ref'};
plots('_C')     = { m.c_id, 'c / mol m^{-3}'};
plots('_logC')  = { logc , 'log(c/mol m^{-3})'};
plots('_absC')  = { absc, '|grad c| / mol m^{-4}'};
plots('_Cx')    = { cx,  'c_x / mol m^{-4}'};
plots('_Cy')    = { cy, 'c_y / mol m^{-4}'};

plots('absn')   = { absn, '|n|'};
plots('nx')     = { m.nx_id,'nx'};
plots('ny')     = { m.ny_id,'ny'};
plots('_absN')   = { absN, '|N|'};
plots('_Nx')     = { m.Nx_id,'Nx'};
plots('_Ny')     = { m.Ny_id,'Ny'};

plots('absDiffusiveFlux') = { absDiffusiveFlux, '-D*|grad c|' };
plots('diffusiveFluxX') = { diffusiveFluxX, '-D*cx' };
plots('diffusiveFluxY') = { diffusiveFluxY, '-D*cy' };

plots('absElectrophoreticFlux') = { absElectrophoreticFlux, '|-z*u*c*F*grad phi|' };
plots('electrophoreticFluxX') = { electrophoreticFluxX, '-z*u*c*F*phix' };
plots('electrophoreticFluxY') = { electrophoreticFluxY, '-z*u*c*F*phiy' };

plots('absTotalFlux') = { absTotalFlux, '|-D*grad c - z*u*c*F*grad phi|' };
plots('totalFluxX') = { totalFluxX, '-D*cx - z*u*c*F*phix' };
plots('totalFluxY') = { totalFluxY, '-D*cy - z*u*c*F*phiy' };


plots('absi')   = { {'sqrt(ix^2+iy^2)'},'|i|'};
plots('ix')     = { {'ix'},'ix'};
plots('iy')     = { {'iy'},'iy'};
% plots('absI')   = { {'sqrt(Ix^2+iy^2)'},'|I|'};
% plots('Ix')     = { {'Ix'},'Ix'};
% plots('Iy')     = { {'Iy'},'Iy'};
plots('_absI')   = { {'tcdee.IlMag'},'|I|'};
plots('_Ix')     = { {'tcdee.Ilx'},'Ix'};
plots('_Iy')     = { {'tcdee.Ily'},'Iy'};

plots('kappa')  = { {'kappa'},'kappa'};

sweepHorizontalCrossection(m,dset,m.L/4,'XleftBoundary','XrightBoundary',0,m.L,plots);
% m.sweepVerticalCrossection(dset,m.W/4);

%% 2d plots
plots = containers.Map;
% plots(title) = { {expression1, expression2, ...}, ylabel };
plots('phi')    = { 'phi/UT', 'phi / U_T'};
plots('gradPhi')   = { {'phix/UT*L','phiy/UT*L'}, 'grad phi * L / U_T'};
% plots('phiy')   = { {'phiy/UT*L'}, 'phi_y * L / U_T'};
plots('_PHI')   = { 'phi', 'phi / V' };
plots('_PHIx')  = { {'phix','phiy'}, 'grad phi / V m^{-1}'};
% plots('_PHIy')  = { {'phiy'}, 'phi_y / V m^{-1}'};

for i=1:m.numberOfSpecies
    plots(m.c_id{i})= { C{i}, 'c / c_ref' };
    plots(sprintf('log%s',m.c_id{i}))   = { logC{i}, 'log(c/c_ref)'};
    plots(sprintf('abs%s',m.c_id{i}))   = { absC{i},'|grad(c)| * L / c_ref'}; 
    plots(sprintf('grad%s',m.cx_id{i}))     = { {Cx{i},Cy{i}},'grad c * L / c_ref'};
%     plots(m.cy_id)     = { Cy, 'c_y * L / c_ref'};
    plots(sprintf('_%s',m.C_id{i}))     = { m.c_id{i}, 'c / mol m^{-3}'};
    plots(sprintf('_log%s',m.C_id{i}))  = { logc{i} , 'log(c/mol m^{-3})'};
    plots(sprintf('_abs%s',m.C_id{i}))  = { absc{i}, '|grad c| / mol m^{-4}'};
    plots(sprintf('_grad%s',m.C_id{i}))    = { {cx{i},cy{i}},  'grad c / mol m^{-4}'};
%     plots(sprintf('_%sy',m.C_id))    = { cy{i}, 'c_y / mol m^{-4}'};

    plots(sprintf('abs%s',m.N_dimless_id{i}))   = { absn{i}, '|n|'};
    plots(sprintf('grad%s',m.N_dimless_id{i}))     = { {m.nx_id{i},m.ny_id{i}},'grad n'};
%     plots(m.ny_id)     = { m.ny_id,'ny'};
    plots(sprintf('_abs%s',m.N_id{i}))   = { absN{i}, '|N|'};
    plots(sprintf('_grad%s',m.N_id{i}))  = { {m.Nx_id{i},m.Ny_id{i}},'grad N'};
%     plots('Ny')     = { m.Ny_id,'Ny'};

    plots(sprintf('absDiffusiveFlux_%s',m.speciesNames{i})) = { absDiffusiveFlux{i}, '-D*|grad c|' };
    plots(sprintf('diffusiveFlux_%s',m.speciesNames{i})) = { {diffusiveFluxX{i},diffusiveFluxY{i}}, '-D*grad c' };
%     plots('diffusiveFluxY') = { diffusiveFluxY, '-D*cy' };

    plots(sprintf('absElectrophoreticFlux_%s',m.speciesNames{i})) = { absElectrophoreticFlux{i}, '|-z*u*c*F*grad phi|' };
    plots(sprintf('electrophoreticFlux_%s',m.speciesNames{i})) = { {electrophoreticFluxX{i},electrophoreticFluxY{i}}, '-z*u*c*F*grad phi' };
%     plots('electrophoreticFluxY') = { electrophoreticFluxY, '-z*u*c*F*phiy' };

    plots(sprintf('absTotalFlux_%s',m.speciesNames{i})) = { absTotalFlux{i}, '|-D*grad c - z*u*c*F*grad phi|' };
    plots(sprintf('totalFluxX_%s',m.speciesNames{i})) = { { totalFluxX{i},totalFluxY{i} }, '-D*grad x - z*u*c*F*grad phi' };
%     plots('totalFluxY') = { totalFluxY, '-z*u*c*F*phiy' };
end

plots('absi')   = { 'sqrt(ix^2+iy^2)','|i|'};
plots('gradi')     = { {'ix','iy'},'grad i'};
% plots('iy')     = { {'iy'},'iy'};
plots('absI')   = { 'sqrt(Ix^2+iy^2)','|I|'};
plots('_gradI')     = { {'Ix','Iy'},'grad I'};
% plots('Iy')     = { {'Iy'},'Iy'};
% plots('absI')   = { {'tcdee.IlMag'},'|I|'};
% plots('Ix')     = { {'tcdee.Ilx'},'Ix'};
% plots('Iy')     = { {'tcdee.Ily'},'Iy'};
plots('kappa')  = { 'kappa','kappa'};

plotStandard2d(m,dset,plots);

%% export
model.result.export.create('exportTertiaryCurrentDistributionData', 'Data');
model.result.export('exportTertiaryCurrentDistributionData').set('data', dset);
model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', m.c_id{i}, i);
end

files('exportTertiaryCurrentDistributionDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistributionDataFile.txt'];

model.result.export('exportTertiaryCurrentDistributionData').set('location', 'grid');
model.result.export('exportTertiaryCurrentDistributionData').set('gridx2', 'range( XleftBoundary, W/1000, XrightBoundary)');
model.result.export('exportTertiaryCurrentDistributionData').set('gridy2', '0');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', files('exportTertiaryCurrentDistributionDataFile'));
model.result.export('exportTertiaryCurrentDistributionData').run;

% full
files('exportTertiaryCurrentDistribution2dDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistribution2dDataFile.txt'];

model.result.export('exportTertiaryCurrentDistributionData').set('location', 'fromdataset');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', files('exportTertiaryCurrentDistribution2dDataFile'));
model.result.export('exportTertiaryCurrentDistributionData').run;

jlh.hf.saveMapAsTxt(files,'globalFiles.txt');
return