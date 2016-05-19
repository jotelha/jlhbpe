function obj = update1dMesh(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    %% meshing
    fprintf('  Setting mesh parameters...\n');
     % get mesh refinement from parameters
%     obj.mesh1D();
    
    obj.m.mesh('standardMesh1d').label('standardMesh');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').selection.geom('geom1d');

    a0_ddl = obj.intFirstDebyeLength(1);
    r_ddl = obj.intFirstDebyeLength(2)/obj.intFirstDebyeLength(1);
    N_ddl = ceil(log(1-obj.epsilon/a0_ddl*(1-r_ddl))/log(r_ddl));

    R_ddl = max(r_ddl^N_ddl, obj.intFirstDebyeLength(end)/obj.intFirstDebyeLength(1));

    a0_eddl = obj.intExtendedDdl(1);
    r_eddl = obj.intExtendedDdl(2)/obj.intExtendedDdl(1);
    N_eddl = ceil(log(1-obj.extendedDdlFactor*obj.epsilon/a0_eddl*(1-r_ddl))/log(r_ddl));

    R_eddl = max(r_eddl^N_eddl, obj.intExtendedDdl(end)/obj.intExtendedDdl(1));
    
    hmax = obj.intRemaining(end);
    
    N_ddl = obj.nvFirstDebyeLength;
    N_eddl = obj.nvExtendedDdl;
    R_ddl = obj.intFirstDebyeLength(end)/obj.intFirstDebyeLength(1);
    R_eddl = obj.intExtendedDdl(end)/obj.intExtendedDdl(1);
    
    
    N_rem = ceil((1-(1+obj.extendedDdlFactor)*obj.epsilon)/obj.intExtendedDdl(end));
    
    %obj.standardMesh.feature('size').set('hauto', 2); % is this the automatic refinement?
    obj.m.mesh('standardMesh1d').feature('size').set('hmax', num2str(hmax) ); % is this the automatic refinement?
    obj.m.mesh('standardMesh1d').feature('size').set('hmin', num2str(a0_ddl) ); 
%     obj.m.mesh('standardMesh1d').feature('size').set('hgrad', num2str(r_zeta) ); 


    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').selection.named('geom1d_ddl1dCumulative_dom');
%     obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('type', 'explicit');
%     obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('explicit', obj.distributionFirstDebyeLengthStr);
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('type','predefined');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('method', 'geometric');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('elemcount', num2str(N_ddl));
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('elemratio', num2str(R_ddl));
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('reverse', 'off');

    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').selection.named('geom1d_extendedDdl1dCumulative_dom');
%     obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('type', 'explicit');
%     obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('explicit', obj.distributionExtendedDdlStr);
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('type','predefined');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('method', 'geometric');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('elemcount', num2str(N_eddl));
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('elemratio', num2str(R_eddl));
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('reverse', 'off');

    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionRemaining1d').selection.named('geom1d_space1dCumulative_dom');
    obj.m.mesh('standardMesh1d').feature('meshingEdge1d').feature('meshDistributionRemaining1d').set('numelem', N_rem);


%     obj.meshDistributionRemaining.selection.named('regionRemaining');
%     obj.meshDistributionRemaining.set('explicit', obj.distributionRemainingStr);
%     obj.meshDistributionRemaining.set('type', 'explicit');

    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').selection.named('geom1d_ddl1dCumulative_dom');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('type','explicit');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionDDL1d').set('explicit', obj.distributionFirstDebyeLengthStr);

    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').selection.named('geom1d_extendedDdl1dCumulative_dom');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('type','explicit');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionExtendedDDL1d').set('explicit', obj.distributionExtendedDdlStr);

    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionRemaining1d').selection.named('geom1d_space1dCumulative_dom');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionRemaining1d').set('type','explicit');
    obj.m.mesh('explicitMesh1d').feature('meshingEdge1d').feature('meshDistributionRemaining1d').set('explicit', obj.distributionRemainingStr);

    obj.m.mesh('standardMesh1d').run;
    obj.m.mesh('explicitMesh1d').run;
end