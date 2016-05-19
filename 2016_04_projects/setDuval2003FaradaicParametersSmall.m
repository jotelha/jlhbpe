% dopeshwarkar2008electrokinetics, small geometry
c = 4;


% concentrations
m.speciesNames = {'H3Op','OHm','Kp','NO3m'};

cKNO3  = [100,10,1,0.1]; % mmol /l = mol/m^3

% % pure aluminium: https://de.wikipedia.org/wiki/Aluminium
% atomicWeightAl = 26.9815 * 1.660539040e-27; %u -> kg
% densityAl = 2.7 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t.
% numberConcentrationAl = densityAl / atomicWeightAl;
% molarConcentrationAl = numberConcentrationAl / constants.AvogadroConstant;
% cAl = molarConcentrationAl;

% pH = 5.5;
pH = 7; % for electroneutrality in bulk
pOH = 14 - pH;
cH3Op  = 10^(-pH)*10^3; %mol/m^3
cOHm   = 10^(-pOH)*10^3;

cKp = cKNO3(c);
cNO3m = cKNO3(c);

% cH2O = 55560; % mol/m^3
% cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html


    
% diffusivities
DKp = 1.89*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Arnold 1955 selfdiffusion
% DKp = 1.957*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
% DNO3m = 1700 / 1e6; % mum^2/s = 10^-6 m^2/s, source http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104439&ver=2
DNO3m = 1.902*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
% DAl3p = 0.541*10^5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
DH3Op = 9.311*1e-5 / 100^2; % source Haynes 2014, 5-77
DOHm = 5.273*1e-5 / 100^2;  % source Haynes 2014, 5-77


% charge numbers
zH3Op = 1;
zOHm = -1;
zKp = 1;
zNO3m = -1;


% switches
m.bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
m.explicitElectrodeGeometry = false;
%bpePoissonEquationBC = 'Dirichlet'; % apply fixed potential at bpe

m.bpeFluxOn = false; % switch species flux due to reactions on or off'

m.widthToHeight = 5/3; % width : height relation for single plots
m.plotsVisible = 'on'; % 'off' stores plots, but does not show them as pop ups


% environment and physical constants
m.epsilon_r = constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

m.T = 298.15; % K, equivalent to 25 deg C

m.numberOfSpecies = 4; 

m.c_bulk  = [cH3Op,cOHm,cKp,cNO3m]; % bulk concentrations
m.z       = [zH3Op,zOHm,zKp,zNO3m]; % charge number of solute species
m.D       = [DH3Op,DOHm,DKp,DNO3m]; % diffusivity of solute species

m.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
m.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD

Wbpe  = 0.5e-3/100; % 500mu m
Wgap  = 0.5e-3/100; %0.5mm / 10000

m.wMeshRatio = 100; % desired ratio between w_mesh and l (l=1)   
m.elementSizeAtSurface = m.epsilon/10;

m.Wbpe = Wbpe;
m.WinsulatorLeft = 0;
m.WinsulatorRight = 0;
m.Wcathode = 0;
m.Wanode = 0;
m.WbulkLeft = Wgap;
m.WbulkRight = Wgap;

m.PHI_bpe     = 0; % [phi] = V
m.PHI_bulk    = 0;

m.deltaPHI = 7.5/100;
m.phiSymmetryFactor = 0.5;

%% reactions

% other paramterized reaction parameters
j0a    = [1.26,3.05,1.85,0.89] / 100; % muA /cm^2 = 100^2 / m^2 * A / 1000^2 = 0.01 A /m^2
j0c    = [0.72,1.36,0.24,0.07] / 100;
% ra = na(1-alpha_a) = n(1-beta), rc = nc alpha_c = n beta
% <=> alpha_a = 1 - ra/na
ra = [33.4,18.4,15.4,10.8]/100;
rc = [28.2,17.7,17.5,14.4]/100;


E_m = [0.095,0.021,-0.133,-0.286];


m.customParameters.cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
m.customParameters.cH2O = 55560; % mol/m^3

m.customParameters.j0a = j0a(c);
m.customParameters.j0c = j0c(c);
m.customParameters.ra = ra(c);
m.customParameters.rc = rc(c);

m.customParameters.E0a = 0.68; % V, duval2003faradaic
m.customParameters.E0c = -0.55;

%% Acidic Oxygen Reduction, balance to anodic oxygen evolution
i = 0;
i = i+1;
m.reactions{i}.name                       = 'AnodicWaterDissotiation';
m.reactions{i}.reaction                   = 'O2 + 4H3O+ + 4e- <= 6H2O'; % simplified H3O+ + e- <=> 3/2 H2O
% could be as well O2 + 4H+ + 4e- <=> 2H2O

% current per reaction
    % Butler-Volmer, cathodic direction > 0 in form of, as in Bard p. 96
    % i = n*F/RT*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
    % in anodic direction > 0:
m.reactions{i}.anodicCurrentTerm        = 'j0a*exp(ra*F_const/RT*(V-E0a))';
m.reactions{i}.cathodicCurrentTerm      = [];
% 
% reactionCurrent = prepTerm(...
%         'n_id*F_const*k0_id*( c_R*exp(n_id*(1-beta_id)*F_const/RT*(V-E0_id)) - c_Ox*exp(-n_id*beta_id*F_const/RT*(V-E0_id)) ) ',...
%         'n_id','k0_id','c_Ox','c_R','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,oxidantConcentrations,reductantConcentrations,obj.beta_id,obj.E0_id);

% m.reactions{i}.oxidants{1}                = 'H3O+'; % consider only water dissotiation
% m.reactions{i}.oxidantConcentrations{1}   = 'cO2*c3'; % because c1 corresponds to cH3Op
% m.reactions{i}.cOx                        = m.customParameters.cO2*cH3Op; % bulk concentrations
% m.reactions{i}.nuOxidants{1}              = -4; 
% m.reactions{i}.oxidantConcentrations{1}   = '0'; 
% m.reactions{i}.cOx                        = 0; % bulk concentrations
% m.reactions{i}.nuOxidants{1}              = 4; 

% m.reactions{i}.oxidants{2}                = 'O2';
% m.reactions{i}.oxidantConcentrations{2} = 'cO2';
% m.reactions{i}.rOxidants{2} = -1;

% m.reactions{i}.reductants{1}              = 'H2O';
% assume only one direction
% m.reactions{i}.reductantConcentrations{1} = 'cH2O';
% m.reactions{i}.reductantConcentrations{1} = '0';
% m.reactions{i}.cR = m.customParameters.cH2O; % assume only one direction, as in Hlushkou 2009
% m.reactions{i}.cR                         = 0; % assume only one direction, as Krawiec
% m.reactions{i}.nuReductants{1}            = 6; % when considering H+, or 6, with h3O+

m.reactions{i}.n                          = 4; % as in Krawiec, but why?
% m.reactions{i}.beta                       = 0.5;
% Zhang, Jiujun 2008 p. 92
% i0 = 2.8e-7 A/cm^2, reference concentrations?
% i0 = 2.8e-7 * 100^2; % A/m^2
% m.reactions{i}.k0 = i0/(F*C0); 
% or, % Pfeiffer p.30 (from 'Redox, Fundamentals, Processes and Applications')
% m.reactions{i}.k0 = 5e-11; % m/s 
% m.reactions{i}.k0                         = 10e-5; % m^4/(s mol), after Krawiec 2008
% m.reactions{i}.E0                         = +1.229; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
m.reactions{i}.flux                       = [-4,0,0,0]; % outward flux of H3Op with cathodic current 

% %% Alkaline oxygen reduction only in anodic reaction, according to Dhopeshwarkar2008electrokinetics
% % or oxygen evolution
% i = i+1;
% m.reactions{i}.name = 'AlkalineOxygenReduction';
% m.reactions{i}.reaction = 'O2 + 2H2O + 4e- <=> 4OH-'; % simplified 1/4 O2 + e- <=> OH-
% 
% m.reactions{i}.oxidants{1} = 'H2O';
% m.reactions{i}.oxidantConcentrations{1} = '0';
% m.reactions{i}.cOx = 0; %only anodic direction
% m.reactions{i}.nuOxidants{1} = -2;
% 
% m.reactions{i}.oxidants{2} = 'O2';
% m.reactions{i}.oxidantConcentrations{2} = '0';
% m.reactions{i}.nuOxidants{2} = -1;
% 
% m.reactions{i}.reductants{1} = 'OH-';
% m.reactions{i}.reductantConcentrations{1} = 'c4'; % because c2 corresponds to cOHm
% m.reactions{i}.cR = cOHm;
% m.reactions{i}.nuReductants{1} = 4;
% 
% m.reactions{i}.n = 4;
% m.reactions{i}.beta = 0.5;
% % Pfeiffer p.30 (from 'Redox, Fundamentals, Processes and Applications')
% m.reactions{i}.k0 = 5e-11; % m/s 
% m.reactions{i}.E0 = 0.401; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
% m.reactions{i}.flux                       = [0,0,0,4,0]; % outward flux of H3Op with cathodic current 


% %% Hydrogen Evolution, only in cathodic direction after  dhopeshwarkar2008electrokinetics
% i = i+1;
% m.reactions{i}.name                   	= 'HydrogenIonReduction'; % assume only happening in one driection
% % m.reactions{i}.reaction = '2H3O+ + 2e- <=> H2 + 2H2O'; % simplified H3O+ + e- <=> 1/2 H2
% m.reactions{i}.reaction                   = '2H+ + 2e- <=> H2'; % only considering H+, not H3O+
% 
% m.reactions{i}.oxidants{1}                = 'H3O+';
% m.reactions{i}.oxidantConcentrations{1}   = 'c1';
% m.reactions{i}.cOx                        = cH3Op; % bulk concentration of hydrogen ions
% m.reactions{i}.nuOxidants{1}              = -2;
% 
% m.reactions{i}.reductants{1}              = 'H2';
% m.reactions{i}.reductantConcentrations{1} = '0';
% % m.reactions{i}.cR                         = m.cH2; % zero
% m.reactions{i}.cR                         = 0; % zero
% m.reactions{i}.nuReductants{1}            = 1;
% 
% % m.reactions{i}.reductants{2} = 'H2O';
% % m.reactions{i}.reductantConcentrations{2} = 'cH2O';
% % m.reactions{i}.nuReductants{2} = 2;
% 
% m.reactions{i}.n                          = 2; % as in Krawiec, why though?
% m.reactions{i}.beta                       = 0.5;
% % Hickling 1940, H2/H+
% % at T = 16 deg C % check again for concentration!
% % data fit to exponentail curve
% % overpotential eta     0.4     0.53    0.64    0.77    V
% % current i             1e-3    1e-2    0.1     1       A/cm^2
% % i0 = 5e-7 A/cm^2
% i0 = 5e-7 * 100^2; % A/m^2
% % m.reactions{i}.k0 = i0 / (2*F*C0);
% m.reactions{i}.k0                         = 1e-10; %m/s % Krawiec 2008
% m.reactions{i}.E0                         = 0; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
% m.reactions{i}.flux                       = [0,0,-2,0,0]; % outward flux of H3Op with cathodic current 

%% after Krawiec 2008, 'Water dissociation', only in cathodic direction
i = i+1;
m.reactions{i}.name                       = 'CathodicWaterDissotiation'; 
% assume only happening in one direction, cathodic
m.reactions{i}.reaction                   = '2H2O + 2e- => H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )

m.reactions{i}.anodicCurrentTerm        = [];
m.reactions{i}.cathodicCurrentTerm      = '-j0c*exp(-rc*F_const/RT*(V-E0c))';

% m.reactions{i}.oxidants{1}                = 'H2O';
% absorbed into rate constant:
% m.reactions{i}.oxidantConcentrations{1} = 'cH2O';
% m.reactions{i}.oxidantConcentrations{1}   = '1';
% m.reactions{i}.cOx = cH2O; 
% m.reactions{i}.cOx                        = 1; % absorbed into rate constant (?)
% m.reactions{i}.nuOxidants{1}              = -2;

% m.reactions{i}.reductants{1}              = 'H2';
% m.reactions{i}.reductantConcentrations{1} = '0';
% % m.reactions{i}.cR                         = m.cH2; % zero
% m.reactions{i}.cR                         = 0; % zero
% m.reactions{i}.nuReductants{1}            = 1;
% m.reactions{i}.reductants{1}              = 'OH-';
% m.reactions{i}.reductantConcentrations{1} = '0';
% m.reactions{i}.cR                         = m.cH2; % zero
% m.reactions{i}.cR                         = 0; % zero
% m.reactions{i}.nuReductants{1}            = 2;
 
% m.reactions{i}.reductants{2} = 'OH-';
% m.reactions{i}.reductantConcentrations{2} = 'c2';
% m.reactions{i}.nuReductants{2} = 2;

m.reactions{i}.n                          = 2; % as in Krawiec
% m.reactions{i}.beta                       = 0.5;
% m.reactions{i}.k0                         = 1e-4; % m^4 / (s mol)
m.reactions{i}.E0                         = -0.8277; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
%m.reactions{i}.flux = [0,1]; % inward flux of OHm with cathodic current 
m.reactions{i}.flux                       = [0,2,0,0]; % inward flux of OHm with cathodic current 

m.nReactions = i;