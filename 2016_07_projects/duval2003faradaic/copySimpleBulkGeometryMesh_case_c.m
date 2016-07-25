% copied from imported mesh part
model.modelNode.create('mcomp1', 'MeshComponent');

model.geom.create('mgeom1', 2); % mesh geometry
model.mesh.create('simpleBulkGeometryRefinedMeshPart', 'mgeom1');
model.mesh('simpleBulkGeometryRefinedMeshPart').create('imp1', 'Import');
model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('source', 'native');
model.mesh('simpleBulkGeometryRefinedMeshPart').feature('imp1').set('filename', 'simpleBulkGeometryRefinedMeshFile');
model.mesh('simpleBulkGeometryRefinedMeshPart').run;

model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
model.mesh('tertiaryCurrentDistributionMesh').create('copy1', 'Copy');

model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').geom(2);
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').geom(2);
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').all;
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_space_dom');
model.mesh('tertiaryCurrentDistributionMesh').run('copy1');