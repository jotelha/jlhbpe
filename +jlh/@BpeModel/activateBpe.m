function obj = setMinimalBC(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
%     fprintf('  Modify boundary conditions of model "%s"...\n',obj.model_tag);
    
    % (de)activate appropriate physics features
    obj.bpeChargeDensityBC.active(true);
    obj.bpePotentialBC.active(false);

%     obj.bulkPotentialBC.active(false);
    
%     obj.BulkConcentrationsBC.active(true); % fixed bulk concentration distribution according do Boltzmann factor
%     obj.FluxAtBulkBoundary.active(false); % no species flux at bulk boundary
    obj.FluxAtSurface.active(true); % no species flux at bulk boundary
%     obj.entireSurfacePotentialGradientBC.active(false);
%     
%     obj.insulatorSurfacePotentialGradientBC.active(false);
%     obj.bulkBoundaryPotentialGradientBC.active(false);
%     
%     obj.zeroSurfaceCurrent.active(false);
  end