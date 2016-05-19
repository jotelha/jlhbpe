function obj = setKrawiecCaseGeneralParameters(obj)
%% set up model compiling and matlab-interna; options

    obj.transformation = 'none'; % apply no transormation to PNP system
    %transformation = 'logConcentration';
    obj.bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
    %bpePoissonEquationBC = 'Dirichlet'; % apply fixed potential at bpe

    obj.bpeFluxOn = false; % switch species flux due to reactions on or off'

    obj.widthToHeight = 5/3; % width : height relation for single plots
    obj.plotsVisible = 'on'; % 'off' stores plots, but does not show them as pop ups


    %% environment and physical constants
    obj.epsilon = 0.01; % relation Debye length : simulation domain, lambdaD / L
    obj.delta = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
    obj.gamma = 0.001; % relation Debye length : simulation domain width, lambdaD / W0


    obj.epsilon_r = constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

    obj.T = 298.15; % K, equivalent to 25 deg C
    
    % pure iron:
%     atomicWeightFe = 55.845 * 1.660539040e-27; %u -> kg, https://en.wikipedia.org/wiki/Iron, https://en.wikipedia.org/wiki/Atomic_mass_unit
%     densityFe = 6.98 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t., https://en.wikipedia.org/wiki/Iron
%     numberConcentrationFe = densityFe / atomicWeightFe;
%     molarConcentrationFe = numberConcentrationFe / constants.AvogadroConstant;

    %% concentrations (bulk or constant)
    % cFe = molarConcentrationFe; % concentration of iron in solid form
    obj.cFe = 1; %absorbed into rate constant
    % obj.cFe2p = 0; % assume constant zero concentration of dissolved ion Fe2+ everywhere
    obj.cH2O = 55560; % mol/m^3
    % saturation concentration of O2 in fresh water at 25 deg C and 1 bar
    obj.cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
    %cO2 = 0.29  % 0.29mM = 0.29e-3 mol/l = 0.29 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html
    obj.cH2 = 0; % assume no hydrogen in solution
    % cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html

    obj.numberOfSpecies = 2;
    pH = 7;
    obj.cH3Op_bulk  = 10^(-pH)*10^3; %mol/m^3
    obj.cOHm_bulk   = obj.cH3Op_bulk;
    obj.cFe2p_bulk  = 0;

    obj.cH3Op   = obj.cH3Op_bulk;
    obj.cOHm    = obj.cOHm_bulk;
    obj.cFe2p   = obj.cFe2p_bulk;
    
    obj.c_bulk  = [obj.cH3Op_bulk,obj.cOHm_bulk]; % bulk concentrations
    obj.z       = [1,-1]; % charge number of solute species
    obj.D       = [jlh.Constants.DiffusivityH3Op,1e-9]; % diffusivity of solute species

    obj.PHI_bpe     = 0.1; % [phi] = V
    obj.PHI_bulk    = 0;
    
    obj.deltaPHI = 0.2;
    obj.phiSymmetryFactor = 0.5;
end