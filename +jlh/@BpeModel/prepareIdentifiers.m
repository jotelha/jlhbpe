%% build identifiers for COMSOL-internal use
function prepareIdentifiers(obj)
    import jlh.hf.*
   
    % dependent on number of species:
%     obj.D_id = nprintf('D%d',obj.numberOfSpecies);
%     obj.z_id = nprintf('z%d',obj.numberOfSpecies);
%     obj.c_id = nprintf('c%d',obj.numberOfSpecies);
%     obj.cx_id = nprintf('c%dx',obj.numberOfSpecies);
%     obj.cy_id = nprintf('c%dy',obj.numberOfSpecies);
%     obj.nx_id = nprintf('n%dx',obj.numberOfSpecies);
%     obj.ny_id = nprintf('n%dy',obj.numberOfSpecies);
%     obj.ix_id = nprintf('i%dx',obj.numberOfSpecies);
%     obj.iy_id = nprintf('i%dy',obj.numberOfSpecies);
% 
%     obj.C_id = nprintf('C%d',obj.numberOfSpecies); % log concentrations
%     obj.c_bulk_id = nprintf('c%d_bulk',obj.numberOfSpecies);
%     obj.c_0_id = nprintf('c%d_0',obj.numberOfSpecies); % initial values throughout domain
%     obj.c_pb_id = nprintf('c%d_pb',obj.numberOfSpecies);
%     obj.cInterpolation_id = nprintf('c%dInterpolation',obj.numberOfSpecies);
%     obj.cSurfaceProbe_id = nprintf('c%dSurfaceProbe',obj.numberOfSpecies);
%     obj.N_id = nprintf('N%d',obj.numberOfSpecies); % species flux at surface
%     obj.NSurfaceProbe_id = nprintf('N%dSurfaceProbe',obj.numberOfSpecies);
%     obj.nu_id_pattern = nprintf('nu%d_%%s',obj.numberOfSpecies); % the stochiometric coefficients of every species for every reaction
% 
%     obj.standardConcentrationsPlot_id = nprintf('standardConcentrationsPlot%d',obj.numberOfSpecies);
%     obj.logConcentrationsPlot_id = nprintf('logConcentrationsPlotGroup%d',obj.numberOfSpecies);

    % 2016-04-26: use species names
    %     obj.speciesNames = speciesNames
    
    obj.D_id = cprintf('D%s',obj.speciesNames);
    obj.z_id = cprintf('z%s',obj.speciesNames);
    obj.c_id = cprintf('c%s',obj.speciesNames);
    obj.cx_id = cprintf('c%sx',obj.speciesNames); % concentration gradient in x direction, dimensionless
    obj.cy_id = cprintf('c%sy',obj.speciesNames);
    obj.nx_id = cprintf('n%sx',obj.speciesNames); % flux in x direction, dimensionless
    obj.ny_id = cprintf('n%sy',obj.speciesNames); % flux in y direction, dimensionless
    obj.ix_id = cprintf('i%sx',obj.speciesNames); % current in x direction, dimensionless
    obj.iy_id = cprintf('i%sy',obj.speciesNames); % current in y direction, dimensionless
    obj.Nx_id = cprintf('N%sx',obj.speciesNames); % flux in x direction 
    obj.Ny_id = cprintf('N%sy',obj.speciesNames); % flux in y direction
    obj.Ix_id = cprintf('I%sx',obj.speciesNames); % current in x direction
    obj.Iy_id = cprintf('I%sy',obj.speciesNames); % current in y direction 

    obj.C_id = cprintf('C%s',obj.speciesNames); % log concentrations
    obj.c_bulk_id = cprintf('c%s_bulk',obj.speciesNames);
    obj.c_0_id = cprintf('c%s_0',obj.speciesNames); % initial values throughout domain
    obj.c_pb_id = cprintf('c%s_pb',obj.speciesNames);
    obj.cInterpolation_id = cprintf('c%sInterpolation',obj.speciesNames);
    obj.cSurfaceProbe_id = cprintf('c%sSurfaceProbe',obj.speciesNames);
    obj.N_id = cprintf('N%s',obj.speciesNames); % species flux at surface
    obj.N_dimless_id = cprintf('NN%s',obj.speciesNames); % species flux at surface, dimensionless
    
    obj.lambdaCBulk_id          = cprintf('lambdaCBulk_%s',obj.speciesNames);
    obj.surfaceFluxBC_id        = cprintf('SurfaceFluxBC_%s',obj.speciesNames);
    obj.bulkFluxBC_id             = cprintf('BulkFluxBC_%s',obj.speciesNames);
    obj.bulkConcentrationBC_id  = cprintf('BulkConcentrationBC_%s',obj.speciesNames);

    obj.NSurfaceProbe_id = cprintf('N%sSurfaceProbe',obj.speciesNames);
    obj.nu_id_pattern = cprintf('nu%s_%%s',obj.speciesNames); % the stochiometric coefficients of every species for every reaction

    obj.standardConcentrationsPlot_id = cprintf('standardConcentrationsPlot%s',obj.speciesNames);
    obj.logConcentrationsPlot_id = cprintf('logConcentrationsPlotGroup%s',obj.speciesNames);

    % dependent on number of surface reactions:
    obj.reactionNames = cellfun( @(r) r.name, obj.reactions', 'UniformOutput', false);
    obj.k0_id = cprintf('k0_%s',obj.reactionNames);
    obj.E0_id = cprintf('E0_%s',obj.reactionNames);
    obj.beta_id = cprintf('beta_%s',obj.reactionNames);
    obj.n_id = cprintf('n_%s',obj.reactionNames);
    obj.i_id = cprintf('i_%s',obj.reactionNames);
    obj.i_dimless_id = cprintf('ii_%s',obj.reactionNames);

    obj.nu_id = cellfun( @(nu) cprintf(nu,obj.reactionNames), obj.nu_id_pattern, 'UniformOutput', false );
    obj.nu_id = [obj.nu_id{:}]';
    obj.iSurfaceProbe_id = cprintf('i%sSurfaceProbe',obj.reactionNames);
    
    % dependent on the bigger one
    obj.multiPurpose1dPlot_id = nprintf('multiPurpose1dPlot%d',obj.nMultiPurpose1dPlots);
end