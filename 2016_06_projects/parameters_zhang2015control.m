% caseStudyTitle = 'duval2003faradaic';
% c = 4;

% concentrations
m.speciesNames = {'H3Op','OHm','DSm','Nap', 'ClO4m'};

% % pure aluminium: https://de.wikipedia.org/wiki/Aluminium
% atomicWeightAl = 26.9815 * 1.660539040e-27; %u -> kg
% densityAl = 2.7 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t.
% numberConcentrationAl = densityAl / atomicWeightAl;
% molarConcentrationAl = numberConcentrationAl /  jlh.Constants.AvogadroConstant;
% cAl = molarConcentrationAl;

% pH = 5.5;
pH = 7; % for electroneutrality in bulk
pOH = 14 - pH;
cH3Op  = 10^(-pH)*10^3; %mol/m^3
cOHm   = 10^(-pOH)*10^3;

cSDS = 1.5e-1; % 0.15mM = 0.15 mol / m^3
cNaClO4 = 1e1; % 0.01 M = 10 mol / m^3

cDSm = cSDS;
cNap = cSDS + cNaClO4;
cClO4m = cNaClO4;

% cH2O = 55560; % mol/m^3
% cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html


    
% diffusivities
DH3Op = 9.311*1e-5 / 100^2; % source Haynes 2014, 5-77
DOHm = 5.273*1e-5 / 100^2;  % source Haynes 2014, 5-77
DDSm = 0.639e-9; % Haynes
DNap = 1.334e-9; % Haynes
DClO4m = 1.792e-9; % Haynes

% charge numbers
zH3Op = 1;
zOHm = -1;
zDSm = -1;
zNap = 1;
zClO4m = -1;

% environment and physical constants
m.epsilon_r =  jlh.Constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

m.T = 298.15; % K, equivalent to 25 deg C

m.numberOfSpecies = 5; 

m.c_bulk  = [cH3Op,cOHm,cDSm,cNap,cClO4m]; % bulk concentrations
m.z       = [zH3Op,zOHm,zDSm,zNap,zClO4m]; % charge number of solute species
m.D       = [DH3Op,DOHm,DDSm,DNap,DClO4m]; % diffusivity of solute species

m.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
m.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD

% Wbpe  = 0.5e-3/100; % 500mu m
% Wgap  = 0.5e-3/100; %0.5mm / 10000

Wbpe  = 1e-2; % 10 mm
Wgap  = 5e-4; % 0.5 mm 
H  = 2e-3; % 2mm, since probing height at 1mm

% m.wMeshRatio = 10; % desired ratio between w_mesh and l (l=1)   
% m.elementSizeAtSurface = m.epsilon/10;

m.Wbpe = Wbpe;
% m.WinsulatorLeft = 0;
% m.WinsulatorRight = 0;
% m.Wcathode = 0;
% m.Wanode = 0;
m.WbulkLeft = Wgap;
m.WbulkRight = Wgap;
m.H = H;

% potentials
% duval2003faradaic kinetic parameters, choose case b (conductivity similar
% to zhang2015control)
c = 2;

E_m = [0.095,0.021,-0.133,-0.286];

m.PHI_bpe     = E_m(c); % [phi] = V
m.PHI_bulk    = 0;

m.deltaPHI = 1; % arbitrary, 1 V
m.phiSymmetryFactor = 0.5;

%% reactions

% other paramterized reaction parameters
j0a    = [1.26,3.05,1.85,0.89] / 100; % muA /cm^2 = 100^2 / m^2 * A / 1000^2 = 0.01 A /m^2
j0c    = [0.72,1.36,0.24,0.07] / 100;
% ra = na(1-alpha_a) = n(1-beta), rc = nc alpha_c = n beta
% <=> alpha_a = 1 - ra/na
ra = [33.4,18.4,15.4,10.8]/100;
rc = [28.2,17.7,17.5,14.4]/100;

% m.customParameters.cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
% m.customParameters.cH2O = 55560; % mol/m^3

m.customParameters.j0a = j0a(c);
m.customParameters.j0c = j0c(c);
m.customParameters.ra = ra(c);
m.customParameters.rc = rc(c);

E0a = 0.68; % V, duval2003faradaic
E0c = -0.55;

% E0a = 1.229; % Haynes 5-83
% E0c = -0.8277; % Haynes

m.customParameters.E0a = E0a;
m.customParameters.E0c = E0c;

%% Acidic Oxygen Reduction, balance to anodic oxygen evolution
i = 0;
i = i+1;
m.reactions{i}.name                       = 'AnodicWaterDissotiation';
m.reactions{i}.reaction                   = 'O2 + 4H3O+ + 4e- <= 6H2O'; % simplified H3O+ + e- <=> 3/2 H2O
% could be as well O2 + 4H+ + 4e- <=> 2H2O

    % in anodic direction > 0:
m.reactions{i}.anodicCurrentTerm        = 'j0a*exp(ra*F_const/RT*(V-E0a))';
m.reactions{i}.cathodicCurrentTerm      = [];

m.reactions{i}.E0                         = E0a;
m.reactions{i}.n                          = 4; 
m.reactions{i}.flux                       = [-4,0,0,0,0]; % outward flux of H3Op with cathodic current 

%% after Krawiec 2008, 'Water dissociation', only in cathodic direction
i = i+1;
m.reactions{i}.name                       = 'CathodicWaterDissotiation'; 
% assume only happening in one direction, cathodic
m.reactions{i}.reaction                   = '2H2O + 2e- => H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )

m.reactions{i}.anodicCurrentTerm        = [];
m.reactions{i}.cathodicCurrentTerm      = '-j0c*exp(-rc*F_const/RT*(V-E0c))';
m.reactions{i}.nuReductants{2} = 2;

m.reactions{i}.E0                         = E0c;
m.reactions{i}.n                          = 2; % as in Krawiec
m.reactions{i}.flux                       = [0,2,0,0,0]; % inward flux of OHm with cathodic current 

m.nReactions = i;

%% simulation settings
smootheningFactor = 100;