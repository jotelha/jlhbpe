function obj = finalizeChoppedMesh(obj)
    obj.m.mesh('standardMesh').create('lateralMesh', 'FreeTri');
    obj.m.mesh('standardMesh').feature('lateralMesh').create('lateralMeshProperties', 'Size');
    obj.m.mesh('standardMesh').feature('lateralMesh').feature('lateralMeshProperties').set('custom', 'on');
    obj.m.mesh('standardMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hmaxactive', 'on');
    obj.m.mesh('standardMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hmax', max(obj.intRemaining));
    obj.m.mesh('standardMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hgradactive', 'on');
    obj.m.mesh('standardMesh').feature('lateralMesh').feature('lateralMeshProperties').set('hgrad', '1.3');
end