function obj = finalizeMappedMesh(obj)
    obj.m.mesh('mappedMesh').create('lateralMesh', 'FreeTri');
    obj.m.mesh('mappedMesh').feature('lateralMesh').create('lateralMeshProperties', 'Size');
    obj.m.mesh('mappedMesh').feature('lateralMesh').feature('lateralMeshProperties').set('custom', 'on');
    obj.m.mesh('mappedMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hmaxactive', 'on');
    obj.m.mesh('mappedMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hmax', max(obj.intRemaining));
    obj.m.mesh('mappedMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hgradactive', 'on');
    obj.m.mesh('mappedMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hgrad', '1.3');
end