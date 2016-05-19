function obj = setUniformFeederElectrodePotentialBC(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    fprintf('  Modify boundary conditions of model "%s"...\n',obj.model_tag);
    
    c0Term =  prepTerm('c_bulk_id/c_ref*exp(-z_id*phi0)','c_bulk_id','z_id',obj.c_bulk_id,obj.z_id); % initial concentration distribution
    phi0Term = '-deltaPhi*x*L/W';
    %% define variables
  
    obj.weVariables.set('phi_s', 'phi_we');
    obj.ceVariables.set('phi_s', 'phi_ce');

    obj.bulkVariables.set('phi_s', phi0Term);
    
    obj.domainVariables.set('phi0', phi0Term);
    for i = 1:obj.numberOfSpecies
        obj.domainVariables.set(obj.c_0_id{i}, c0Term{i});
    end
 
    % (de)activate appropriate physics features
    obj.bpeChargeDensityBC.active(true);
    obj.bulkPotentialBC.active(false);
    obj.bpePotentialBC.active(false);
    
    obj.BulkConcentrationsBC.active(true); % fixed bulk concentration distribution according do Boltzmann factor
    obj.FluxAtBulkBoundary.active(true); % no species flux at bulk boundary
    obj.entireSurfacePotentialGradientBC.active(true);
    obj.bulkBoundaryPotentialGradientBC.active(true);
    
    obj.zeroSurfaceCurrent.active(false);
  end