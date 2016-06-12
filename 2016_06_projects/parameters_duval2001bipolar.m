% caseStudyTitle = 'duval2001bipolar';
c = 4;
% aluminium electrode
% 10^-2 M KNO3 solution
% at phi < -1.6 V, reduction of water sets in
% at phi > 1.3 V, anodic dissolution of aluminium sets in
% for KNO3 c = 2*10^-4 ~ 10^-1 M, current i <= 0.2mA cm^-2
% for KNO3 c = 10^-1 M dissolution current  10 mA cm^-2 at 2.5V
% voltammograms from -1 to -2.5 V and from 2.5 to - 1 V
% E0a = -1.82 V Ag | Ag | KCl (sat) for Al -> Al3+ + 3e- at pH = 5.5
% E0c = -0.55 V Ag | Ag | KCl(sat) for 2H2O + 2e- -> H2 + 2OH-
% reference electrode Ag | Ag | KCl(sat) +0.222 mV vs SHE
% c/ mmol l^-1  j0a/muA cm^-2   ra * 10^2   j0c/muA cm^-2   rc * 10^2   Em / V
% 100           0.08            7           2.12            10.7        -0.58
% 50            0.10            6.5         1.57            10.3        -0.62
% 10            0.13            5.6         1.16            9.3         -0.65
% 0.2           5.32            1.6         0.21            6.9         -1.76
m.speciesNames = {'H3Op','OHm','Kp','NO3m'};

% concentrations
cKNO3  = [100,50,10,0.2]; % mmol /l = mol/m^3

% pure aluminium: https://de.wikipedia.org/wiki/Aluminium
atomicWeightAl = 26.9815 * 1.660539040e-27; %u -> kg
densityAl = 2.7 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t.
numberConcentrationAl = densityAl / atomicWeightAl;
molarConcentrationAl = numberConcentrationAl / jlh.Constants.AvogadroConstant;
cAl = molarConcentrationAl;

pH = 7;
pOH = 14 - pH;
cH3Op  = 10^(-pH)*10^3; %mol/m^3
cOHm   = 10^(-pOH)*10^3;


cKp = cKNO3(c);
cNO3m = cKNO3(c);

cH2O = 55560; % mol/m^3
% cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3,
% http://www1.lsbu.ac.uk/water/electrolysis.html'

% conductivity estimations:
% % KNO3 in aqueous solution
% xiKNO3_haynes   = [0.5, 1, 2, 5, 10, 15, 20, 25];
% % kappa in mS/cm concentration in mass percent
% %   c       0.5     1       2       5       10      15      20      25
% %   kappa   5.5     10.7    20.1    47.0    87.3    124     157     182
% % source: Electrical conductivity of aqueous solutions
% %   Hayne 2014 CRC Handbook of chemistry and physics, 5-73
% 
% cKNO3_haynes    = [0.05, 0.099, 0.2, 0.509, 1.051, (1.509+1.747)/2, 2.24, 2.759] * 10^3; % mol / L = 10^3 mol/m^3
% % source: Electrical conductivity of aqueous solutions
% %   Hayne 2014 CRC Handbook of chemistry and physics, 5-73
% kappa_haynes = [5.5,10.7,20.1,47,87.3,124,157,182] *0.1; %mS/cm = 100/1000 S/m
% 
% % molar conductivity
% cKNO3_haynes_dilute = [0.0005,0.001,0.005,0.01,0.02,0.05,0.1]*10^3; % mol/L = 10^3 mol/m^3
% Lambda_haynes_dilute = [142.7,141.77,138.41,132.75,132.34,126.25,120.34] * 1e-4; %m^2 S / mol
% % source: equiValent conductiVity of electrolytes in aqueous solution
% %   Hayne 2014 CRC Handbook of chemistry and physics, 5-76
% kappa_haynes_dilute = Lambda_haynes_dilute.*cKNO3_haynes_dilute; % Lambda = kappa / c
% 
% % cKNO3_interp = [cKNO3_haynes_dilute,cKNO3_haynes(3:end)];
% % kappa_interp = [kappa_haynes_dilute,kappa_haynes(3:end)];
% 
% % kappa = interp1(cKNO3_haynes_dilute,kappa_haynes_dilute,cKNO3);
% % extrapolation by quadratic function:
% coeffs = polyfit(cKNO3_haynes_dilute,kappa_haynes_dilute,2);
% kappa = polyval(coeffs,cKNO3);


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

% measures
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
E_m = [-0.58,-0.62,-0.65,-1.76];

m.PHI_bpe     = E_m(c); % [phi] = V
m.PHI_bulk    = 0;

m.deltaPHI = 7.5/100;
m.phiSymmetryFactor = 0.5;

%% reactions

m.customParameters.cAl   = cAl;
m.customParameters.cH2O  = cH2O;

% paramterized reaction parameters
j0a    = [0.08,0.10,0.13,5.32] / 100; % muA /cm^2 = 100^2 / m^2 * A / 1000^2 = 0.01 A /m^2
j0c    = [2.12,1.57,1.16,0.21] / 100;
% ra = na(1-alpha_a) = n(1-beta), rc = nc alpha_c = n beta
% <=> alpha_a = 1 - ra/na
ra = [7,6.5,5.6,1.6]/100;
rc = [10.7,10.3,9.3,6.9]/100;

E0a = -1.82; % V, duval2001bipolar
E0c = -0.55;

m.customParameters.j0a = j0a(c);
m.customParameters.j0c = j0c(c);
m.customParameters.ra = ra(c);
m.customParameters.rc = rc(c);

m.customParameters.E0a = E0a;
m.customParameters.E0c = E0c;

%% Anodic Aluminium Dissolution
i = 0;
i = i+1;
m.reactions{i}.name                       = 'AnodicAluminiumDissolution';
m.reactions{i}.reaction                   = 'Al3+ + 3e- <= Al'; 

% current per reaction
    % Butler-Volmer, cathodic direction > 0 in form of, as in Bard p. 96
    % i = n*F/RT*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
    % in anodic direction > 0:
m.reactions{i}.anodicCurrentTerm        = 'j0a*exp(ra*F_const/RT*(V-E0a))';
m.reactions{i}.cathodicCurrentTerm      = [];

m.reactions{i}.n                          = 3;
m.reactions{i}.E0                         = E0a; % V, source: duval2001bipolar
m.reactions{i}.flux                       = [0,0,0,0]; % no flux

%% after Krawiec 2008, 'Water dissociation', only in cathodic direction
i = i+1;
m.reactions{i}.name                       = 'CathodicWaterDissotiation'; 
% assume only happening in one direction, cathodic
m.reactions{i}.reaction                   = '2H2O + 2e- => H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )

m.reactions{i}.anodicCurrentTerm        = [];
m.reactions{i}.cathodicCurrentTerm      = '-j0c*exp(-rc*F_const/RT*(V-E0c))';

m.reactions{i}.n                          = 2; 
m.reactions{i}.E0                         = E0c; % V, source: duval2001bipolar
m.reactions{i}.flux                       = [0,2,0,0]; % inward flux of OHm with cathodic current 

m.nReactions = i;