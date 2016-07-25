import com.comsol.model.*
import com.comsol.model.util.*
import jlh.*
import jlh.hf.*

if ~exist('caseStudyParameterFile','var')
    caseStudyParameterFile = 'parameters_duval2001bipolar.m';
end

caseStudyTitle = 'pde1d';
createEmptyProject;

% or, to load an existing project from server or from file
% tag = 'xyz'
% loadExistingProject


%% get file names
% T = readtable('geometryPartsFile.txt','ReadRowNames',true,'ReadVariableNames',false,'Delimiter',' ');
T = readtable('globalFiles.txt','ReadRowNames',false,'ReadVariableNames',false,'Delimiter',' ');
% writetable(T,'geometryPartsFile.txt','WriteRowNames',true,'WriteVariableNames',false,'Delimiter',' ');
C = table2cell(T);
files = containers.Map(C(:,1),C(:,2));
% geometryPartsMphFile = T.(1)('geometryPartsMphFile');
% simpleBulkGeometryRefinedMeshFile = T.(1)('simpleBulkGeometryRefinedMeshFile');
% simpleDdlGeometryRefinedMeshFile = T.(1)('simpleDdlGeometryRefinedMeshFile');
%% load parameters, geometry and mesh parts
model = m.m;

model.param().loadFile(files('parameterFile'));

% mesh parts

model.modelNode.create('mcomp1', 'MeshComponent');
model.geom.create('mgeom1', 1); % mesh geometry
model.mesh.create('simple1dGeometryRefinedMeshPart', 'mgeom1');
model.mesh('simple1dGeometryRefinedMeshPart').create('imp1', 'Import');
model.mesh('simple1dGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
model.mesh('simple1dGeometryRefinedMeshPart').feature('imp1').set('filename', files('simple1dGeometryRefinedMeshFile'));
model.mesh('simple1dGeometryRefinedMeshPart').run;

%%
% makeDilutedSpeciesAndElectrostatics1dComponent
makeCoefficientPDE1dComponent

pname = {'X'};
plistarr = {'range(XleftBoundary,W/100,XrightBoundary)'};
makeSimpleStudy
% tags of created features in sol_id, study_id, studyStep_id, compile_id,
% variables_id, solver_id and store_id

% adaption for 1d study

model.study(study_id).feature(studyStep_id).set('pname', pname);
model.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
model.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname); 
model.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);

model.study(study_id).feature(studyStep_id).set('preusesol', 'yes');

model.sol(sol_id).feature(solver_id).set('stol', '1e-14');
model.sol(sol_id).feature(solver_id).set('keeplog', 'on');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'auto');
model.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');

model.sol(sol_id).run(store_id);

storedSolutions = model.study(study_id).getSolverSequences('Stored');
solInfo = mphsolinfo(model,'soltag',char(storedSolutions(1)));
dset = char(solInfo.dataset);

% make parametric extrusion
model.result.dataset.create('par1', 'Parametric1D');
model.result.dataset('par1').set('data', dset);
% model.result.dataset('par1').setIndex('looplevelinput', 'manual', 2);
% model.result.dataset('par1').setIndex('looplevelinput', 'first', 1);
% model.result.dataset('par1').setIndex('looplevelinput', 'manual', 1);
model.result.dataset('par1').set('levelscaleactive', 'on');
model.result.dataset('par1').set('levelscale', '1');
model.result.dataset('par1').set('innerinput', 'all');
model.result.dataset('par1').run;

m.saveState;

createPlots

%%
% set up expressions to plot
% small letters: dimensionles, capital letters: dimensional

logc = prepTerm('log(c_id)','c_id',m.c_id);
cx = prepTerm('c_idx','c_id',m.c_id);
% cy = prepTerm('tds.grad_c_idy','c_id',m.c_id);
% absc  = prepTerm('sqrt(cx^2+cy^2)','cx','cy',cx,cy);
% dimensional
C = prepTerm('c_id/c_ref','c_id',m.c_id);
logC = prepTerm('log(c_id/c_ref)','c_id',m.c_id);
Cx = prepTerm('c_idx/c_ref*L','c_id',m.c_id);
% Cy = prepTerm('(tds.grad_c_idy/c_ref*L)','c_id',m.c_id);
% absC  = prepTerm('sqrt(Cx^2+Cy^2)','Cx','Cy',Cx,Cy);

% absn = prepTerm('sqrt(nx_id^2+ny_id^2)','nx_id','ny_id',m.nx_id,m.ny_id);
% absN = prepTerm('sqrt(Nx_id^2+Ny_id^2)','Nx_id','Ny_id',m.Nx_id,m.Ny_id);

% absDiffusiveFlux = prepTerm('tds.dfluxMag_c_id','c_id',m.c_id);
diffusiveFluxX = prepTerm('-D_id*c_idx','D_id','c_id',m.D_id,m.c_id);
% diffusiveFluxY = prepTerm('tds.dflux_c_idy','c_id',m.c_id);
electrophoreticFluxX = prepTerm('-z_id*D_id*c_id*F_const/RT*phix','z_id','D_id','c_id',m.z_id,m.D_id,m.c_id);
% electrophoreticFluxY = prepTerm('tds.mflux_c_idy','c_id',m.c_id);
% absElectrophoreticFlux = prepTerm('sqrt(efx^2+efy^2)','efx','efy',electrophoreticFluxX,electrophoreticFluxY);
% absTotalFlux = prepTerm('tds.tfluxMag_c_id','c_id',m.c_id);
totalFluxX = prepTerm('diffusiveFluxX+electrophoreticFluxX','diffusiveFluxX','electrophoreticFluxX',diffusiveFluxX,electrophoreticFluxX);
% totalFluxY = prepTerm('tds.tflux_c_idy','c_id',m.c_id);
       
            
plots = containers.Map;
% plots(title) = { {expression1, expression2, ...}, ylabel };
plots('phi')    = { {'phi/UT'}, 'phi / U_T'};
plots('phix')   = { {'phix/UT*L'}, 'phi_x * L / U_T'};
% plots('phiy')   = { {'phiy/UT*L'}, 'phi_y * L / U_T'};
plots('_PHI')   = { {'phi'}, 'phi / V' };
plots('_PHIx')  = { {'phix'}, 'phi_x / V m^{-1}'};
% plots('_PHIy')  = { {'phiy'}, 'phi_y / V m^{-1}'};

plots('c')      = { C, 'c / c_ref' };
plots('logc')   = { logC, 'log(c/c_ref)'};
% plots('absc')   = { absC,'|grad(c)| * L / c_ref'}; 
plots('cx')     = { Cx,'c_x * L / c_ref'};
% plots('cy')     = { Cy, 'c_y * L / c_ref'};
plots('_C')     = { m.c_id, 'c / mol m^{-3}'};
plots('_logC')  = { logc , 'log(c/mol m^{-3})'};
% plots('_absC')  = { absc, '|grad c| / mol m^{-4}'};
plots('_Cx')    = { cx,  'c_x / mol m^{-4}'};
% plots('_Cy')    = { cy, 'c_y / mol m^{-4}'};

% plots('absn')   = { absn, '|n|'};
plots('nx')     = { m.nx_id,'nx'};
% plots('ny')     = { m.ny_id,'ny'};
% plots('absN')   = { absN, '|N|'};
plots('Nx')     = { m.Nx_id,'Nx'};
% plots('Ny')     = { m.Ny_id,'Ny'};

% plots('absDiffusiveFlux') = { absDiffusiveFlux, '-D*|grad c|' };
plots('diffusiveFluxX') = { diffusiveFluxX, '-D*cx' };
% plots('diffusiveFluxY') = { diffusiveFluxY, '-D*cy' };

% plots('absElectrophoreticFlux') = { absElectrophoreticFlux, '|-z*u*c*F*grad phi|' };
plots('electrophoreticFluxX') = { electrophoreticFluxX, '-z*u*c*F*phix' };
% plots('electrophoreticFluxY') = { electrophoreticFluxY, '-z*u*c*F*phiy' };

% plots('absTotalFlux') = { absTotalFlux, '|-D*grad c - z*u*c*F*grad phi|' };
plots('totalFluxX') = { totalFluxX, '-D*cx - z*u*c*F*phix' };
% plots('totalFluxY') = { totalFluxY, '-z*u*c*F*phiy' };


% plots('absi')   = { {'sqrt(ix^2+iy^2)'},'|i|'};
plots('ix')     = { {'ix'},'ix'};
% plots('iy')     = { {'iy'},'iy'};
% plots('absI')   = { {'sqrt(Ix^2+iy^2)'},'|I|'};
plots('Ix')     = { {'Ix'},'Ix'};
% plots('Iy')     = { {'Iy'},'Iy'};
% plots('absI')   = { {'tcdee.IlMag'},'|I|'};
% plots('Ix')     = { {'tcdee.Ilx'},'Ix'};
% plots('Iy')     = { {'tcdee.Ily'},'Iy'};

plots('kappa')  = { {'kappa'},'kappa'};

% XleftBoundary = model.param.evaluate('XleftBoundary');
% XrightBoundary = model.param.evaluate('XrightBoundary');
% sweepHorizontalCrossection(m,dset,m.W/4,'XleftBoundary','XrightBoundary',XleftBoundary,XrightBoundary,plots);
% m.sweepVerticalCrossection(dset,m.W/4);

plotParametric1d(m,dset,dset,'all',plots);


%% surface plots, global
plots = containers.Map;

c = prepTerm('intSurface(c_id)','c_id',m.c_id);
logc = prepTerm('intSurface(log(c_id))','c_id',m.c_id);
cx = prepTerm('intSurface(c_idx)','c_id',m.c_id);
C = prepTerm('intSurface(c_id/c_ref)','c_id',m.c_id);
logC = prepTerm('intSurface(log(c_id/c_ref))','c_id',m.c_id);
Cx = prepTerm('intSurface(c_idx/c_ref*L)','c_id',m.c_id);
diffusiveFluxX = prepTerm('intSurface(diffusiveFluxX)','diffusiveFluxX',diffusiveFluxX);
electrophoreticFluxX = prepTerm('intSurface(electrophoreticFluxX)','electrophoreticFluxX',electrophoreticFluxX);
totalFluxX = prepTerm('intSurface(totalFluxX)','totalFluxX',totalFluxX);
nx = prepTerm('intSurface(nx_id)','nx_id',m.nx_id);
Nx = prepTerm('intSurface(Nx_id)','Nx_id',m.Nx_id);

plots('phi')    = { {'intSurface(phi/UT)'}, 'phi / U_T'};
plots('phix')   = { {'intSurface(phix/UT*L)'}, 'phi_x * L / U_T'};
plots('_PHI')   = { {'intSurface(phi)'}, 'phi / V' };
plots('_PHIx')  = { {'intSurface(phix)'}, 'phi_x / V m^{-1}'};

plots('c')      = { C, 'c / c_ref' };
plots('logc')   = { logC, 'log(c/c_ref)'};
plots('cx')     = { Cx,'c_x * L / c_ref'};
plots('_C')     = { c, 'c / mol m^{-3}'};
plots('_logC')  = { logc , 'log(c/mol m^{-3})'};
plots('_Cx')    = { cx,  'c_x / mol m^{-4}'};

plots('nx')     = { nx,'nx'};
plots('Nx')     = { Nx,'Nx'};
plots('diffusiveFluxX') = { diffusiveFluxX, '-D*cx' };
plots('electrophoreticFluxX') = { electrophoreticFluxX, '-z*u*c*F*phix' };
plots('totalFluxX') = { totalFluxX, '-D*cx - z*u*c*F*phix' };

plots('ix')     = { {'intSurface(ix)'},'ix'};
plots('Ix')     = { {'intSurface(Ix)'},'Ix'};


plots('kappa')  = { {'intSurface(kappa)'},'kappa'};

plots('i_total')    = { {'intSurface(i_total)'}, 'I / A m^-2'};
plots('i_cathodic')   = { {'intSurface(i_cathodic)'}, 'I / A m^-2'};
plots('i_anodic')   = { {'intSurface(i_anodic)'}, 'I / A m^-2' };
plots('log_i_cathodic')   = { {'log(abs(intSurface(i_cathodic)))'}, 'I / A m^-2'};
plots('log_i_anodic')   = { {'log(abs(intSurface(i_anodic)))'}, 'I / A m^-2' };
plots('log_i')   = { {'log(abs(intSurface(i_cathodic)))','log(abs(intSurface(i_anodic)))'}, 'I / A m^-2' };
for i = 1:m.nReactions
        plots(m.i_id{i}) = { { sprintf('intSurface(%s)',m.i_id{i}) },  'I / A m^-2'};
end
for i = 1:m.numberOfSpecies
    plots(m.N_id{i}) = { { sprintf('intSurface(%s)',m.N_id{i}) },  'mol / s m^-2'};
    plots(m.Nx_id{i})= { { sprintf('intSurface(%s)',m.Nx_id{i}) }, 'mol / s m^-2'};
end

for i = 1:m.numberOfSpecies
    plots([m.N_id{i},'_comparative']) = { { sprintf('intSurface(%s)',m.N_id{i}), sprintf('intSurface(%s)',m.Nx_id{i}),},  'mol / s m^-2'};
end

plotGlobal1d(m,dset,'X',plots);
%% export
model.result.export.create('exportDilutedSpeciesAndElectrostatics1dData', 'Data');
model.result.export('exportDilutedSpeciesAndElectrostatics1dData').set('data', 'par1');
model.result.export('exportDilutedSpeciesAndElectrostatics1dData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportDilutedSpeciesAndElectrostatics1dData').setIndex('expr', m.c_id{i}, i);
end

files('exportDilutedSpeciesAndElectrostatics1dDataFile') = [pwd,'\',m.projectPath,'\exportDilutedSpeciesAndElectrostatics1dDataFile.txt'];

model.result.export('exportDilutedSpeciesAndElectrostatics1dData').set('location', 'fromdataset');
% model.result.export('exportDilutedSpeciesAndElectrostatics1dData').set('gridx2', 'range( XleftBoundary, W/1000, XrightBoundary)');
% model.result.export('exportDilutedSpeciesAndElectrostatics1dData').set('gridy2', '0');
model.result.export('exportDilutedSpeciesAndElectrostatics1dData').set('filename', files('exportDilutedSpeciesAndElectrostatics1dDataFile'));
model.result.export('exportDilutedSpeciesAndElectrostatics1dData').run;

jlh.hf.saveMapAsTxt(files,'globalFiles.txt');
return