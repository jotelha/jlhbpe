function createProbes(obj)
%  obj.surfaceProbeTable = obj.m.result.table.create('probeTable', 'Table');
% 
%     obj.phiSurfaceProbe = obj.m.result.numerical.create('phiSurfaceProbe', 'EvalPoint');
% 
%     obj.cSurfaceProbe = cell(obj.numberOfSpecies,1);
%     obj.NSurfaceProbe = cell(obj.numberOfSpecies,1);
%     for i=1:obj.numberOfSpecies  
%         obj.cSurfaceProbe{i} = obj.m.result.numerical.create(obj.cSurfaceProbe_id{i}, 'EvalPoint');
%         obj.NSurfaceProbe{i} = obj.m.result.numerical.create(obj.NSurfaceProbe_id{i}, 'EvalPoint');
%     end
% 
%     obj.iCathodicSurfaceProbe = obj.m.result.numerical.create('iCathodic', 'EvalPoint');
%     obj.iAnodicSurfaceProbe = obj.m.result.numerical.create('iAnodic', 'EvalPoint');
%     obj.iTotalSurfaceProbe = obj.m.result.numerical.create('iTotal', 'EvalPoint');
    obj.surfaceProbeTable = obj.m.result.table.create('probeTable', 'Table');

%     obj.phiSurfaceProbe = obj.m.result.numerical.create('phiSurfaceProbe', 'EvalGlobal');

%     obj.cSurfaceProbe = cell(obj.numberOfSpecies,1);
    obj.NSurfaceProbe = cell(obj.numberOfSpecies,1);
    for i=1:obj.numberOfSpecies  
%         obj.cSurfaceProbe{i} = obj.m.result.numerical.create(obj.cSurfaceProbe_id{i}, 'EvalGlobal');
        obj.NSurfaceProbe{i} = obj.m.result.numerical.create(obj.NSurfaceProbe_id{i}, 'EvalGlobal');
    end
    for i=1:obj.nReactions
        obj.iSurfaceProbe{i} = obj.m.result.numerical.create(obj.iSurfaceProbe_id{i}, 'EvalGlobal');
    end
    
    

    obj.iCathodicSurfaceProbe = obj.m.result.numerical.create('iCathodic', 'EvalGlobal');
    obj.iAnodicSurfaceProbe = obj.m.result.numerical.create('iAnodic', 'EvalGlobal');
    obj.iTotalSurfaceProbe = obj.m.result.numerical.create('iTotal', 'EvalGlobal');
end