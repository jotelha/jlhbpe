% a simplemost mesh
model.mesh.create('tertiaryCurrentDistributionMesh', 'tertiaryCurrentDistributionGeometry');
model.mesh('tertiaryCurrentDistributionMesh').create('ftri1', 'FreeTri');

model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('custom', 'on');
model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('hgrad', '1.2162996682304414');
model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('hmin', '1.751100e-03*L');
model.mesh('tertiaryCurrentDistributionMesh').feature('size').set('hmax', 'H');
model.mesh('tertiaryCurrentDistributionMesh').feature('ftri1').set('yscale', '20');
model.mesh('tertiaryCurrentDistributionMesh').run;