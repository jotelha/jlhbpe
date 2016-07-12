% copied from imported mesh part
model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
model.mesh('tertiaryCurrentDistributionMesh').create('copy1', 'Copy');

model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').set('mesh', 'simpleBulkGeometryRefinedMeshPart');
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').geom(2);
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').geom(2);
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('source').all;
model.mesh('tertiaryCurrentDistributionMesh').feature('copy1').selection('destination').named('tertiaryCurrentDistributionGeometry_simpleBulkGeometryPartInstance1_space_dom');
model.mesh('tertiaryCurrentDistributionMesh').run('copy1');