function obj = replicateTestingMeshPrototype(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    for i=2:obj.nMeshChops
        id = sprintf('replicateMeshPrototype%d',i-1);
        dest1 = 2*i;
        dest2 = 2*i+1;
        obj.m.mesh('manualMesh').create(id, 'Copy');
    %     model.mesh('manualMesh').feature('copy1').selection('destination').named('geom_meshChopPrototype_dom');
        obj.m.mesh('manualMesh').feature(id).selection('source').named('geom_unitCellPrototypeCumulative_dom');
%         obj.m.mesh('manualMesh').feature(id).selection('source').set([2, 3});
        obj.m.mesh('manualMesh').feature(id).selection('destination').set([dest1, dest2]);
    end
end

