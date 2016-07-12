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

if ~exist('BatchDir','var')
    BatchDir = [ pwd, '\', m.projectPath, '\batch'];
end

if ~exist(BatchDir,'dir')
    mkdir(BatchDir);
end

GlobalEvaluationsTableFileName = [BatchDir,'\globalEvaluations.txt'];

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
meshFile = 'tertiaryCurrentDistributionMesh_copy.m'; %% standard method of meshing, copy mesh from geometry file
% meshFile = 'tertiaryCurrentDistributionMesh_simple.m'; %% tells following script to run this file for meshing

makeTertiaryCurrentDistribution2dComponent

% makeBatchStudy
% tags of created features in sol_id, study_id, studyStep_id, compile_id,
% variables_id, solver_id and store_id
% model.sol('sol1').feature('s1').feature('dDef').set('ooc', 'on');

%% numerical evaluations
model.result.table.create('tbl1', 'Table');
% model.result.table('tbl1').comments('Global Evaluation 1 (intWE(tcdee.Ilx))');
model.result.export.create('tbl1', 'Table');
model.result.export('tbl1').set('filename', GlobalEvaluationsTableFileName);


titles = {'PHI_bpe', 'Ix_WE', 'Ix_CE', 'Iy_BPE', ' I_total', 'I_anodic', 'I_cathodic', 'I_faradaic', 'I_ohmic'};
expressions = { 'PHI_bpe', ... % mixed potential
                'intWE(tcdee.Ilx)', ...
                'intCE(tcdee.Ilx)', ...
                'intBPE(tcdee.Ily)', ...
                'intBPE(i_total)',...
                'intBPE(i_anodic)',...
                'intBPE(i_cathodic)', ...
                '(intBPE(i_anodic)-intBPE(i_cathodic))/2', ... % faradaic current, averaged
                '(intWE(tcdee.Ilx)+intCE(tcdee.Ilx))/2 - (intBPE(i_anodic)-intBPE(i_cathodic))/2' }; % ohmic current, averaged
            
for i=1:numel(expressions)
    model.result.numerical.create(titles{i}, 'EvalGlobal');
    model.result.numerical(titles{i}).set('probetag', 'none');
    model.result.numerical(titles{i}).set('table', 'tbl1');
    model.result.numerical(titles{i}).set('expr', expressions{i});
    model.result.numerical(titles{i}).set('descr', expressions{i});
end



%% export
model.result.export.create('exportTertiaryCurrentDistributionData', 'Data');
model.result.export('exportTertiaryCurrentDistributionData').set('data', dset);
model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', 'phi', 0);
for i=1:m.numberOfSpecies
    model.result.export('exportTertiaryCurrentDistributionData').setIndex('expr', m.c_id{i}, i);
end

% files('exportTertiaryCurrentDistributionDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistributionDataFile.txt'];
files('exportTertiaryCurrentDistributionDataFile') = [BatchDir,'\exportTertiaryCurrentDistributionDataFile.txt'];

model.result.export('exportTertiaryCurrentDistributionData').set('location', 'grid');
model.result.export('exportTertiaryCurrentDistributionData').set('gridx2', 'range( XleftBoundary, W/1000, XrightBoundary)');
model.result.export('exportTertiaryCurrentDistributionData').set('gridy2', '0');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', files('exportTertiaryCurrentDistributionDataFile'));
model.result.export('exportTertiaryCurrentDistributionData').run;

% full
% files('exportTertiaryCurrentDistribution2dDataFile') = [pwd,'\',m.projectPath,'\exportTertiaryCurrentDistribution2dDataFile.txt'];
files('exportTertiaryCurrentDistribution2dDataFile') = [BatchDir,'\exportTertiaryCurrentDistribution2dDataFile.txt'];

model.result.export('exportTertiaryCurrentDistributionData').set('location', 'fromdataset');
model.result.export('exportTertiaryCurrentDistributionData').set('filename', files('exportTertiaryCurrentDistribution2dDataFile'));
% model.result.export('exportTertiaryCurrentDistributionData').run;

%% make batch study
pname = { 'DeltaPHI' };
plistarr = { '1.5,2.5,3.5,4.5' };

makeBatchStudy

m.saveState;