function obj = replicateMeshPrototype(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    for i=2:obj.nMeshChops
        id = sprintf('replicateMeshPrototype%d',i-1);
        dest1 = 2*i;
        dest2 = 2*i+1;
        obj.m.mesh('standardMesh').create(id, 'Copy');
    %     model.mesh('standardMesh').feature('copy1').selection('destination').named('geom_meshChopPrototype_dom');
        obj.m.mesh('standardMesh').feature(id).selection('source').named('geom_unitCellPrototypeCumulative_dom');
%         obj.m.mesh('standardMesh').feature(id).selection('source').set([2, 3});
        obj.m.mesh('standardMesh').feature(id).selection('destination').set([dest1, dest2]);
    end
end

