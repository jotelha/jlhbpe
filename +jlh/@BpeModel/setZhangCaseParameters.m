function obj = setZhangCaseParameters(obj)
    obj.bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
    %bpePoissonEquationBC = 'Dirichlet'; % apply fixed potential at bpe

    obj.bpeFluxOn = false; % switch species flux due to reactions on or off'

    obj.widthToHeight = 5/3; % width : height relation for single plots
    obj.plotsVisible = 'on'; % 'off' stores plots, but does not show them as pop ups


    %% environment and physical constants
%     obj.epsilon = 0.01; % relation Debye length : simulation domain, lambdaD / L
%     obj.delta = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
%     obj.gamma = 0.001; % relation Debye length : simulation domain width, lambdaD / W0


    obj.epsilon_r = constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

    obj.T = 298.15; % K, equivalent to 25 deg C
    
    % pure iron:
%     atomicWeightFe = 55.845 * 1.660539040e-27; %u -> kg, https://en.wikipedia.org/wiki/Iron, https://en.wikipedia.org/wiki/Atomic_mass_unit
%     densityFe = 6.98 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t., https://en.wikipedia.org/wiki/Iron
%     numberConcentrationFe = densityFe / atomicWeightFe;
%     molarConcentrationFe = numberConcentrationFe / constants.AvogadroConstant;

    %% concentrations (bulk or constant)
    % cFe = molarConcentrationFe; % concentration of iron in solid form
    cFe = 1; %absorbed into rate constant
    % obj.cFe2p = 0; % assume constant zero concentration of dissolved ion Fe2+ everywhere
    cH2O = 55560; % mol/m^3
    % saturation concentration of O2 in fresh water at 25 deg C and 1 bar
    cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
    %cO2 = 0.29  % 0.29mM = 0.29e-3 mol/l = 0.29 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html
    % cH2 = 0; % assume no hydrogen in solution
    cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html

    obj.numberOfSpecies = 5;
    pH = 7;
    cH3Op_bulk  = 10^(-pH)*10^3; %mol/m^3
    cOHm_bulk   = obj.cH3Op_bulk;
    cFe2p_bulk  = 0;

    cH3Op   = cH3Op_bulk;
    cOHm    = cOHm_bulk;
    cFe2p   = cFe2p_bulk;
    
    % concentration of SDS: 0.15 mM
    cSDS_bulk = 0.15; % mM = 1/1000 mol/l = M/m^3
    % concentration of NaClO4: 0.01 M
    cNaClO4_bulk = 0.01 * 10^3; % 1 M = 1 mol/l = 1000 mol/m^3
    
    cNap_bulk   = cSDS_bulk+cNaClO4_bulk;
    cDSm_bulk   = cSDS_bulk;
    cClO4m_bulk = cNaClO4_bulk;

    % H3O+,OH-,Na+,DS-,ClO4-
    obj.c_bulk  = [cH3Op_bulk,cOHm_bulk,cNap_bulk,cDSm_bulk,cClO4m_bulk]; % bulk concentrations
    obj.z       = [1,-1,1,-1,-1]; % charge number of solute species
    obj.D       = [jlh.Constants.DiffusivityH3Op,1e-9,1e-9,1e-9,1e-9]; % diffusivity of solute species

    Wbpe  = 0.01; % 10mm
    Wgap  = 0.0005; %0.5mm
    W     = Wbpe+2*Wgap;
    L     = Wbpe/10;
    
    obj.epsilon     = obj.lambdaD/L; % relation Debye length : simulation domain, lambdaD / L
    obj.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
    obj.gamma       = obj.lambdaD/W; % relation Debye length : simulation domain width, lambdaD / W0
    obj.beta        = obj.lambdaD/Wbpe;

    obj.PHI_bpe     = 0.1; % [phi] = V
    obj.PHI_bulk    = 0;
    
    obj.deltaPHI = 0.2;
    obj.phiSymmetryFactor = 0.5;
    
      %% redox reactions Ox + n e- <=> R
    % k0 = 1e-3; % standard rate constant, 1 cm/s = 1e-3 m/s

    % chemists need to approve concentration boundary conditions
    %% iron ion reduction
    i = 1;                
%     obj.reactions{i}.name                       = 'IronIonReduction';
%     obj.reactions{i}.reaction                   = 'Fe2+ + 2e- <=> Fe';
%     obj.reactions{i}.oxidants{1}                = 'Fe2+';
%     obj.reactions{i}.oxidantConcentrations{1}   = 'cFe2p';
%     obj.reactions{i}.cOx                        = obj.cFe2p; % zero
%     obj.reactions{i}.nuOxidants{1}              = -1;
%     obj.reactions{i}.reductants{1}              = 'Fe';
%     obj.reactions{i}.reductantConcentrations{1} = 'cFe';
%     % obj.reactions{i}.cR = cFe;
%     obj.reactions{i}.cR                         = 1; % absorbed into rate constant
%     obj.reactions{i}.nuReductants{1}            = 1;
%     obj.reactions{i}.n                          = 2;
%     obj.reactions{i}.beta                       = 0.5;
%     % Tanaka 1964 p.975, Fe/Fe2+ 
%     % for T = 25 deg C and COx = CR = C = 1M,  
%     % log i0 = 4.046, i0 in A/cm^2
%     i0 = 10^4.046 * 100^2; % A/m^2
%     C0 = 10^3; % 1 M = 1 mol/l = 10^3 mol/m^3
%     % obj.reactions{i}.k0 = i0 / (2*F*C0);
%     obj.reactions{i}.k0                         = 3e-10; % mol / (s m^2), after Krawiec 2008
%     obj.reactions{i}.E0                         = -0.44; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
%     obj.reactions{i}.flux                       = [0,0]; % no flux of H3Op and OHm, should be in the format of -z

    %% Acidic Oxygen Reduction
%     i = i+1;
    obj.reactions{i}.name                       = 'AcidicOxygenReduction';
    obj.reactions{i}.reaction                   = 'O2 + 4H3O+ + 4e- <=> 6H2O'; % simplified H3O+ + e- <=> 3/2 H2O
    % could be as well O2 + 4H+ + 4e- <=> 2H2O
    obj.reactions{i}.oxidants{1}                = 'H3O+';
    obj.reactions{i}.oxidantConcentrations{1}   = 'cO2*c1'; % because c1 corresponds to cH3Op
    obj.reactions{i}.cOx                        = cO2*cH3Op; % bulk concentrations
    obj.reactions{i}.nuOxidants{1}              = -4; 

    obj.reactions{i}.oxidants{2}                = 'O2';
    % obj.reactions{i}.oxidantConcentrations{2} = 'cO2';
    % obj.reactions{i}.rOxidants{2} = -1;

    obj.reactions{i}.reductants{1}              = 'H2O';
    obj.reactions{i}.reductantConcentrations{1} = 'cH2O';
%     obj.reactions{i}.reductantConcentrations{1} = '0';
    obj.reactions{i}.cR = cH2O;
%     obj.reactions{i}.cR                         = 0; % assume only one direction, as Krawiec
    obj.reactions{i}.nuReductants{1}            = 2; % when considering H+, or 6, with h3O+

    obj.reactions{i}.n                          = 4; % as in Krawiec, but why?
    obj.reactions{i}.beta                       = 0.5;
    % Zhang, Jiujun 2008 p. 92
    % i0 = 2.8e-7 A/cm^2, reference concentrations?
    % i0 = 2.8e-7 * 100^2; % A/m^2
    % obj.reactions{i}.k0 = i0/(F*C0); 
    % or, % Pfeiffer p.30 (from 'Redox, Fundamentals, Processes and Applications')
    % obj.reactions{i}.k0 = 5e-11; % m/s 
    obj.reactions{i}.k0                         = 10e-5; % m^4/(s mol), after Krawiec 2008
    obj.reactions{i}.E0                         = +1.229; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    obj.reactions{i}.flux                       = [-4,0,0,0,0]; % outward flux of H3Op with cathodic current 

    %% Alkaline oxygen reduction
    i = i+1;
    obj.reactions{i}.name = 'AlkalineOxygenReduction';
    obj.reactions{i}.reaction = 'O2 + 2H2O + 4e- <=> 4OH-'; % simplified 1/4 O2 + e- <=> OH-
  
    obj.reactions{i}.oxidants{2} = 'O2';
    obj.reactions{i}.oxidantConcentrations{2} = 'cO2';
    obj.reactions{i}.nuOxidants{2} = -1;
  
    obj.reactions{i}.oxidants{1} = 'H2O';
    obj.reactions{i}.oxidantConcentrations{1} = 'cH2O';
    obj.reactions{i}.nuOxidants{2} = -2;
    obj.reactions{i}.cOx = cH2O*cO2;
    
    obj.reactions{i}.reductants{1} = 'OH-';
    obj.reactions{i}.reductantConcentrations{1} = 'c2'; % because c2 corresponds to cOHm
    obj.reactions{i}.cR = cOHm_bulk;
    obj.reactions{i}.nuReductants{1} = 4;
    obj.reactions{i}.n = 4;
    obj.reactions{i}.beta = 0.5;
    % Pfeiffer p.30 (from 'Redox, Fundamentals, Processes and Applications')
    obj.reactions{i}.k0 = 5e-11; % m/s 
    obj.reactions{i}.E0 = -0.401; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    obj.reactions{i}.flux                       = [0,4,0,0,0]; % outward flux of H3Op with cathodic current 

    
    %% Hydrogen Evolution
    i = i+1;
    obj.reactions{i}.name                   	= 'HydrogenIonReduction'; % assume only happening in one driection
    % obj.reactions{i}.reaction = '2H3O+ + 2e- <=> H2 + 2H2O'; % simplified H3O+ + e- <=> 1/2 H2
    obj.reactions{i}.reaction                   = '2H+ + 2e- <=> H2'; % only considering H+, not H3O+

    obj.reactions{i}.oxidants{1}                = 'H3O+';
    obj.reactions{i}.oxidantConcentrations{1}   = 'c1';
    obj.reactions{i}.cOx                        = cH3Op; % bulk concentration of hydrogen ions
    obj.reactions{i}.nuOxidants{1}              = -2;

    obj.reactions{i}.reductants{1}              = 'H2';
    obj.reactions{i}.reductantConcentrations{1} = 'cH2';
    obj.reactions{i}.cR                         = cH2; % zero
    obj.reactions{i}.nuReductants{1}            = 1;

    % obj.reactions{i}.reductants{2} = 'H2O';
    % obj.reactions{i}.reductantConcentrations{2} = 'cH2O';
    % obj.reactions{i}.nuReductants{2} = 2;

    obj.reactions{i}.n                          = 2; % as in Krawiec, why though?
    obj.reactions{i}.beta                       = 0.5;
    % Hickling 1940, H2/H+
    % at T = 16 deg C % check again for concentration!
    % data fit to exponentail curve
    % overpotential eta     0.4     0.53    0.64    0.77    V
    % current i             1e-3    1e-2    0.1     1       A/cm^2
    % i0 = 5e-7 A/cm^2
    i0 = 5e-7 * 100^2; % A/m^2
    % obj.reactions{i}.k0 = i0 / (2*F*C0);
    obj.reactions{i}.k0                         = 1e-10; %m/s % Krawiec 2008
    obj.reactions{i}.E0                         = 0; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    obj.reactions{i}.flux                       = [-2,0,0,0,0]; % outward flux of H3Op with cathodic current 

    %% after Krawiec 2008, 'Water dissociation'
    i = i+1;
    obj.reactions{i}.name                       = 'WaterReduction'; 
    % assume only happening in one direction, cathodic
    obj.reactions{i}.reaction                   = '2H2O + 2e- <=> H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )

    obj.reactions{i}.oxidants{1}                = 'H2O';
    % absorbed into rate constant:
    %obj.reactions{i}.oxidantConcentrations{1} = 'cH2O';
    obj.reactions{i}.oxidantConcentrations{1}   = '1';
    % obj.reactions{i}.cOx = cH2O; 
    obj.reactions{i}.cOx                        = 1; % absorbed into rate constant (?)
    obj.reactions{i}.nuOxidants{1}              = -2;

    obj.reactions{i}.reductants{1}              = 'H2';
    obj.reactions{i}.reductantConcentrations{1} = 'cH2';
    obj.reactions{i}.cR                         = cH2; % zero
    obj.reactions{i}.nuReductants{1}            = 1;

    % obj.reactions{i}.reductants{2} = 'OH-';
    % obj.reactions{i}.reductantConcentrations{2} = 'c2';
    % obj.reactions{i}.nuReductants{2} = 2;

    obj.reactions{i}.n                          = 2; % as in Krawiec, not sure abouyt it
    obj.reactions{i}.beta                       = 0.5;
    obj.reactions{i}.k0                         = 1e-4; % m^4 / (s mol)
    obj.reactions{i}.E0                         = -0.8277; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    %obj.reactions{i}.flux = [0,1]; % inward flux of OHm with cathodic current 
    obj.reactions{i}.flux                       = [0,2,0,0,0]; % inward flux of OHm with cathodic current 

    obj.nReactions = i;
end