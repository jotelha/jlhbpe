function obj = update1dFluxOperators(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% local definitions
    fprintf('  Setting up local definitions...\n');
   
    % integration WE
    obj.integrateWE.selection.named('workingElectrode');
    obj.integrateWE.label('integrateWE');
    obj.integrateWE.set('opname', 'integrateWE');

    % integrate CE
    obj.integrateCE.selection.named('counterElectrode');
    obj.integrateCE.label('integrateCE');
    obj.integrateCE.set('opname', 'integrateCE');
    
    % integrate domain
    obj.integrateDomain.selection.all;
    obj.integrateDomain.label('integrateDomain');
    obj.integrateDomain.set('opname', 'integrateDomain');

    % maximum on whole computational domain
    obj.maximumOnDomain.set('opname', 'maximumOnDomain');
    obj.maximumOnDomain.selection.all;
end