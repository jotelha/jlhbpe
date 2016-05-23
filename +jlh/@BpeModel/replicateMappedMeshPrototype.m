function obj = replicateMappedMeshPrototype(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
%     for i=2:obj.nMeshChops
    if obj.nMeshChops >= 4
        prototypeDomains = double(obj.m.selection('geom_unitCellPrototypeCumulative_dom').entities(2));        
        for i=3:(obj.nMeshChops-1)
            id = sprintf('replicateMappedMeshPrototype%d',i-1);
%             id_upper = sprintf('replicateMappedMeshPrototype%d_upper',i-1);
%             id_lower = sprintf('replicateMappedMeshPrototype%d_lower',i-1);

            destinationDomains = prototypeDomains + (i-2)*[2;2;1];
%             dest1 = 3*i+1;
%             dest2 = 3*i+2;
%             dest3 = 3*i+3;
            obj.m.mesh('mappedMesh').create(id, 'Copy');
%             model.mesh('mappedMesh').feature('copy1').selection('destination').named('geom_meshChopPrototype_dom');
            obj.m.mesh('mappedMesh').feature(id).selection('source').named('geom_unitCellPrototypeCumulative_dom');
            obj.m.mesh('mappedMesh').feature(id).selection('destination').set(destinationDomains);
    
%             obj.m.mesh('mappedMesh').create(id_lower, 'CopyDomain');
%             obj.m.mesh('mappedMesh').create(id_upper, 'CopyDomain');
%             obj.m.mesh('mappedMesh').feature(id_lower).selection('source').set([7 8]);
%             obj.m.mesh('mappedMesh').feature(id_upper).selection('source').set(9);
%             obj.m.mesh('mappedMesh').feature(id_lower).selection('destination').set([dest1, dest2]);
%             obj.m.mesh('mappedMesh').feature(id_upper).selection('destination').set(dest3);
        end
    end
end

