function createEvaluations1d(obj)
%     obj.surfaceProbeTable = obj.m.result.table.create('probeTable', 'Table');
% 
%     obj.NSurfaceProbe = cell(obj.numberOfSpecies,1);
%     for i=1:obj.numberOfSpecies  
% %         obj.cSurfaceProbe{i} = obj.m.result.numerical.create(obj.cSurfaceProbe_id{i}, 'EvalGlobal');
%         obj.NSurfaceProbe{i} = obj.m.result.numerical.create(obj.NSurfaceProbe_id{i}, 'EvalGlobal');
%     end
%     for i=1:obj.nReactions
%         obj.iSurfaceProbe{i} = obj.m.result.numerical.create(obj.iSurfaceProbe_id{i}, 'EvalGlobal');
%     end
%     
%     
% 
%     obj.iCathodicSurfaceProbe = obj.m.result.numerical.create('iCathodic', 'EvalGlobal');
%     obj.iAnodicSurfaceProbe = obj.m.result.numerical.create('iAnodic', 'EvalGlobal');
%     obj.iTotalSurfaceProbe = obj.m.result.numerical.create('iTotal', 'EvalGlobal');
    import jlh.*;
    import jlh.hf.*;
    
    obj.m.result.table.create('multiPurposeEvaluationTable1d', 'Table');

    obj.m.result.table.create('surfaceEvaluationTable1d', 'Table');
    obj.m.result.table.create('bulkExitEvaluationTable1d', 'Table');


    de = flattenCell(obj.domainExpressions1d);
    dl = obj.e2t(de);
%     dl = strrep(de,'(','_');
%     dl = strrep(dl,')','');
%     dl = strrep(dl,'/','_by_');
%     dl = strrep(dl,'*','_times_');
%     dl = strrep(dl,'-','_minus_');
%     dl = strrep(dl,'+','_plus_');


    se = flattenCell(obj.surfaceExpressions1d);
    sl = obj.e2t(se);
%     sl = strrep(de,'(','_');
%     sl = strrep(sl,')','');
%     sl = strrep(sl,'/','_by_');
%     sl = strrep(sl,'*','_times_');
%     sl = strrep(sl,'-','_minus_');
%     sl = strrep(sl,'+','_plus_');
    
    
    for i = 1:numel(de)     
        lSurface = sprintf('%s_surface',dl{i});
        obj.m.result.numerical.create(lSurface, 'EvalPoint');
        obj.m.result.numerical(lSurface).label(lSurface);
        obj.m.result.numerical(lSurface).set('expr',de{i});
        obj.m.result.numerical(lSurface).set('table','surfaceEvaluationTable1d');
%         obj.m.result.numerical(lSurface).selection.named('surfaceVertex1d');
        
        lBulk = sprintf('%s_bulk',dl{i});
        obj.m.result.numerical.create(lBulk, 'EvalPoint');
        obj.m.result.numerical(lBulk).label(lBulk);
        obj.m.result.numerical(lBulk).set('expr',de{i});
        obj.m.result.numerical(lBulk).set('table','bulkExitEvaluationTable1d');
%         obj.m.result.numerical(lBulk).selection.named('surfaceVertex1d');
    end
    for i = 1:numel(se)
        obj.m.result.numerical.create(sl{i}, 'EvalPoint');
        obj.m.result.numerical(sl{i}).label(sl{i});
        obj.m.result.numerical(sl{i}).set('expr',se{i});
        obj.m.result.numerical(sl{i}).set('table','surfaceEvaluationTable1d');
%         obj.m.result.numerical(sl{i}).selection.named('surfaceVertex1d');
    end
        
end