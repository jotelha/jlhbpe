% a simplemost mesh
model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
model.mesh('tertiaryCurrentDistributionMesh').create('ftri1', 'FreeTri');

model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('custom', 'on');
model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('hgrad', '1.1705007841642914');
model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('hmin', '1.284265e-03*L');
model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('hmax', 'H');
model.mesh('tertiaryCurrentDistributionMesh').feature('ftri1').set('yscale', '20');
model.mesh('tertiaryCurrentDistributionMesh').run;