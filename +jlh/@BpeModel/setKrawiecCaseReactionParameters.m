function obj = setKrawiecCaseReactionParameters(obj)
    %% redox reactions Ox + n e- <=> R
    % k0 = 1e-3; % standard rate constant, 1 cm/s = 1e-3 m/s

    % chemists need to approve concentration boundary conditions
    %% iron ion reduction
    i = 1;                
    obj.reactions{i}.name                       = 'IronIonReduction';
    obj.reactions{i}.reaction                   = 'Fe2+ + 2e- <=> Fe';
    obj.reactions{i}.oxidants{1}                = 'Fe2+';
    obj.reactions{i}.oxidantConcentrations{1}   = 'cFe2p';
    obj.reactions{i}.cOx                        = obj.cFe2p; % zero
    obj.reactions{i}.nuOxidants{1}              = -1;
    obj.reactions{i}.reductants{1}              = 'Fe';
    obj.reactions{i}.reductantConcentrations{1} = 'cFe';
    % obj.reactions{i}.cR = cFe;
    obj.reactions{i}.cR                         = 1; % absorbed into rate constant
    obj.reactions{i}.nuReductants{1}            = 1;
    obj.reactions{i}.n                          = 2;
    obj.reactions{i}.beta                       = 0.5;
    % Tanaka 1964 p.975, Fe/Fe2+ 
    % for T = 25 deg C and COx = CR = C = 1M,  
    % log i0 = 4.046, i0 in A/cm^2
    i0 = 10^4.046 * 100^2; % A/m^2
    C0 = 10^3; % 1 M = 1 mol/l = 10^3 mol/m^3
    % obj.reactions{i}.k0 = i0 / (2*F*C0);
    obj.reactions{i}.k0                         = 3e-10; % mol / (s m^2), after Krawiec 2008
    obj.reactions{i}.E0                         = -0.44; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    obj.reactions{i}.flux                       = [0,0]; % no flux of H3Op and OHm, should be in the format of -z

    %% Acidic Oxygen Reduction
    i = i+1;
    obj.reactions{i}.name                       = 'AcidicOxygenReduction';
    obj.reactions{i}.reaction                   = 'O2 + 4H3O+ + 4e- <=> 6H2O'; % simplified H3O+ + e- <=> 3/2 H2O
    % could be as well O2 + 4H+ + 4e- <=> 2H2O
    obj.reactions{i}.oxidants{1}                = 'H3O+';
    obj.reactions{i}.oxidantConcentrations{1}   = 'cO2*c1'; % because c1 corresponds to cH3Op
    obj.reactions{i}.cOx                        = obj.cO2*obj.cH3Op; % bulk concentrations
    obj.reactions{i}.nuOxidants{1}              = -4; 

    obj.reactions{i}.oxidants{2}                = 'O2';
    % obj.reactions{i}.oxidantConcentrations{2} = 'cO2';
    % obj.reactions{i}.rOxidants{2} = -1;

    obj.reactions{i}.reductants{1}              = 'H2O';
    % assume only one direction
    % obj.reactions{i}.reductantConcentrations{1} = 'cH2O';
    obj.reactions{i}.reductantConcentrations{1} = '0';
    % obj.reactions{i}.cR = cH2O;
    obj.reactions{i}.cR                         = 0; % assume only one direction, as Krawiec
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
    obj.reactions{i}.flux                       = [-4,0]; % outward flux of H3Op with cathodic current 

    %% Alkaline oxygen reduction
    % i = i+1;
    % obj.reactions{i}.name = 'AlkalineOxygenReduction';
    % obj.reactions{i}.reaction = 'O2 + 2H2O + 4e- <=> 4OH-'; % simplified 1/4 O2 + e- <=> OH-
    % obj.reactions{i}.oxidants{2} = 'O2';
    % obj.reactions{i}.oxidantConcentrations{2} = 'cO2';
    % obj.reactions{i}.rOxidants{2} = 1/4;
    % obj.reactions{i}.oxidants{1} = 'H2O';
    % obj.reactions{i}.oxidantConcentrations{1} = 'cH2O';
    % obj.reactions{i}.cOx = cH2O;
    % obj.reactions{i}.nuOxidants{1} = 1/2;
    % obj.reactions{i}.reductants{1} = 'OH-';
    % obj.reactions{i}.reductantConcentrations{1} = 'c2'; % because c2 corresponds to cOHm
    % obj.reactions{i}.cR = cOHm_bulk;
    % obj.reactions{i}.nuReductants{1} = 1;
    % obj.reactions{i}.n = 1;
    % obj.reactions{i}.beta = 0.5;
    % % Pfeiffer p.30 (from 'Redox, Fundamentals, Processes and Applications')
    % obj.reactions{i}.k0 = 5e-11; % m/s 
    % obj.reactions{i}.E0 = -0.401; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)

    %% Hydrogen Evolution
    i = i+1;
    obj.reactions{i}.name                   	= 'HydrogenIonReduction'; % assume only happening in one driection
    % obj.reactions{i}.reaction = '2H3O+ + 2e- <=> H2 + 2H2O'; % simplified H3O+ + e- <=> 1/2 H2
    obj.reactions{i}.reaction                   = '2H+ + 2e- <=> H2'; % only considering H+, not H3O+

    obj.reactions{i}.oxidants{1}                = 'H3O+';
    obj.reactions{i}.oxidantConcentrations{1}   = 'c1';
    obj.reactions{i}.cOx                        = obj.cH3Op; % bulk concentration of hydrogen ions
    obj.reactions{i}.nuOxidants{1}              = -2;

    obj.reactions{i}.reductants{1}              = 'H2';
    obj.reactions{i}.reductantConcentrations{1} = 'cH2';
    obj.reactions{i}.cR                         = obj.cH2; % zero
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
    obj.reactions{i}.flux                       = [-2,0]; % outward flux of H3Op with cathodic current 

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
    obj.reactions{i}.cR                         = obj.cH2; % zero
    obj.reactions{i}.nuReductants{1}            = 1;

    % obj.reactions{i}.reductants{2} = 'OH-';
    % obj.reactions{i}.reductantConcentrations{2} = 'c2';
    % obj.reactions{i}.nuReductants{2} = 2;

    obj.reactions{i}.n                          = 2; % as in Krawiec, not sure abouyt it
    obj.reactions{i}.beta                       = 0.5;
    obj.reactions{i}.k0                         = 1e-4; % m^4 / (s mol)
    obj.reactions{i}.E0                         = -0.8277; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    %obj.reactions{i}.flux = [0,1]; % inward flux of OHm with cathodic current 
    obj.reactions{i}.flux                       = [0,0]; % inward flux of OHm with cathodic current 

    obj.nReactions = i;
end