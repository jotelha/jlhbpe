function updateEvaluations1d(obj,dset,eval,table_id)
    import jlh.*;
    import jlh.hf.*;
    
    % eval is cell array of probes to be evaluated
    
%     obj.m.result.table.create('surfaceEvaluationTable1d', 'Table');
%     obj.m.result.table.create('bulkExitEvaluationTable1d', 'Table');

    if( ~exist('table_id','var') )
        table_id = 'multiPurposeEvaluationTable1d';
    end

    obj.m.result.table(table_id).clearTableData();
%     obj.m.result.table('surfaceEvaluationTable1d').clearTableData();
%     obj.m.result.table('bulkExitEvaluationTable1d').clearTableData();
    
    de = flattenCell(obj.domainExpressions1d);
    dl = obj.e2t(de);


    se = flattenCell(obj.surfaceExpressions1d);
    sl = obj.e2t(se);
    
    if( ~exist('eval','var') )
        eval = 'all';
    end

    for i = 1:numel(de)    
        lSurface = sprintf('%s_surface',dl{i});
        if( any(strcmp('all',eval)) || any(strcmp(dl{i},eval)))
            obj.m.result.numerical(lSurface).set('data',dset);
            obj.m.result.numerical(lSurface).selection.named('surfaceVertex1d');
    %         obj.m.result.numerical(lSurface).setResult;
            obj.m.result.numerical(lSurface).set('table',table_id);
            if(obj.m.result.table(table_id).getNRows == 0)
                            obj.m.result.numerical(lSurface).setResult;
            else
                obj.m.result.numerical(lSurface).appendResult;
            end
        end

        lBulk = sprintf('%s_bulk',dl{i});
        if( any(strcmp(eval,'all')) || any(strcmp(dl{i},eval)))
            obj.m.result.numerical(lBulk).set('data',dset);
            obj.m.result.numerical(lBulk).selection.named('bulkExitVertex1d');
            obj.m.result.numerical(lBulk).set('table',table_id);
            if(obj.m.result.table(table_id).getNRows == 0)
                obj.m.result.numerical(lBulk).setResult;
            else
                obj.m.result.numerical(lBulk).appendResult;
            end
        end
    end
    for i = 1:numel(se)
        if( any(strcmp(eval,'all')) || any(strcmp(sl{i},eval)))
            obj.m.result.numerical(sl{i}).set('data',dset);
            obj.m.result.numerical(sl{i}).selection.named('surfaceVertex1d');
            obj.m.result.numerical(sl{i}).set('table',table_id);
            if(obj.m.result.table(table_id).getNRows == 0)
                obj.m.result.numerical(sl{i}).setResult;
            else
                obj.m.result.numerical(sl{i}).appendResult;
            end
        end
    end
        
end