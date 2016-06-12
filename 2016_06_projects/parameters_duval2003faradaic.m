% caseStudyTitle = 'duval2003faradaic';
c = 4;

% concentrations
m.speciesNames = {'H3Op','OHm','Kp','NO3m'};

cKNO3  = [100,10,1,0.1]; % mmol /l = mol/m^3

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
m.epsilon_r =  jlh.Constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

m.T = 298.15; % K, equivalent to 25 deg C

m.numberOfSpecies = 4; 

m.c_bulk  = [cH3Op,cOHm,cKp,cNO3m]; % bulk concentrations
m.z       = [zH3Op,zOHm,zKp,zNO3m]; % charge number of solute species
m.D       = [DH3Op,DOHm,DKp,DNO3m]; % diffusivity of solute species

m.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
m.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD

Wbpe  = 0.5e-3/100; % 500mu m
Wgap  = 0.5e-3/100; %0.5mm / 10000

m.wMeshRatio = 10; % desired ratio between w_mesh and l (l=1)   
m.elementSizeAtSurface = m.epsilon/10;

m.Wbpe = Wbpe;
m.WinsulatorLeft = 0;
m.WinsulatorRight = 0;
m.Wcathode = 0;
m.Wanode = 0;
m.WbulkLeft = Wgap;
m.WbulkRight = Wgap;

% potentials
E_m = [0.095,0.021,-0.133,-0.286];

m.PHI_bpe     = E_m(c); % [phi] = V
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

m.customParameters.cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
m.customParameters.cH2O = 55560; % mol/m^3

m.customParameters.j0a = j0a(c);
m.customParameters.j0c = j0c(c);
m.customParameters.ra = ra(c);
m.customParameters.rc = rc(c);

E0a = 0.68; % V, duval2003faradaic
E0c = -0.55;
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
m.reactions{i}.n                          = 4; % as in Krawiec, but why?
m.reactions{i}.flux                       = [-4,0,0,0]; % outward flux of H3Op with cathodic current 

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
m.reactions{i}.flux                       = [0,2,0,0]; % inward flux of OHm with cathodic current 

m.nReactions = i;