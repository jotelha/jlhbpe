% caseStudyTitle = 'duval2003faradaic';

% concentrations
m.speciesNames = {'Al3p','H3Op','OHm','Kp','NO3m'};

kappa = 10; % 10^-1 Ohm^-1 cm^-1 = 0.1 S /cm * 100 cm / m = 10 S / m

% relate conductivity with concentrations
%% electrical conductivity:
% KNO3 in aqueous solution
% kappa in mS/cm concentration in mass percent
%   c       0.5     1       2       5       10      15      20      25
%   kappa   5.5     10.7    20.1    47.0    87.3    124     157     182
% source: Electrical conductivity of aqueous solutions
%   Hayne 2014 CRC Handbook of chemistry and physics, 5-73
xiKNO3_kappa_haynes   = [0.5, 1, 2, 5, 10, 15, 20, 25]; % in mass percent
kappa_haynes = [5.5,10.7,20.1,47,87.3,124,157,182] *0.1; % in S/m = 1 / (Ohm m) %mS/cm = 100/1000 S/m

xiKNO3_c_haynes = [0.5,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24];
cKNO3_haynes    = [0.05, 0.099, 0.2, 0.302, 0.405, 0.509, 0.615, 0.722, 0.830, 0.940, 1.051, 1.277, 1.509, 1.747, 1.991, 2.240, 2.497, 2.759] * 10^3; % mol / L = 10^3 mol/m^3
% source: Electrical conductivity of aqueous solutions
%   Hayne 2014 CRC Handbook of chemistry and physics, 5-138

xiKNO3_from_c = @(c) interp1(cKNO3_haynes,xiKNO3_c_haynes,c);
kappaKNO3_from_xi = @(xi) interp1(xiKNO3_kappa_haynes,kappa_haynes,xi);
kappaKNO3 = @(c) kappaKNO3_from_xi( xiKNO3_from_c(c) );

xiKNO3_from_kappa = @(kappa) interp1(kappa_haynes,xiKNO3_kappa_haynes,kappa);
cKNO3_from_xi = @(xi) interp1(xiKNO3_c_haynes,cKNO3_haynes,xi);
cKNO3_from_kappa = @(kappa) cKNO3_from_xi( xiKNO3_from_kappa(kappa) );

% molar conductivity
cKNO3_haynes_dilute = [0.0005,0.001,0.005,0.01,0.02,0.05,0.1]*10^3; % mol/L = 10^3 mol/m^3
Lambda_haynes_dilute = [142.7,141.77,138.41,132.75,132.34,126.25,120.34] * 1e-4; %m^2 S / mol
% source: equiValent conductiVity of electrolytes in aqueous solution
%   Hayne 2014 CRC Handbook of chemistry and physics, 5-76
kappa_haynes_dilute = Lambda_haynes_dilute.*cKNO3_haynes_dilute; % Lambda = kappa / c

% cKNO3_interp = [cKNO3_haynes_dilute,cKNO3_haynes(3:end)];
% kappa_interp = [kappa_haynes_dilute,kappa_haynes(3:end)];

% kappa = interp1(cKNO3_haynes_dilute,kappa_haynes_dilute,cKNO3);
% extrapolation by quadratic function:
% coeffs = polyfit(cKNO3_haynes_dilute,kappa_haynes_dilute,2);
% kappa = polyval(coeffs,cKNO3);


% cKNO3 = [ 7.19e2 ] ;

% % pure aluminium: https://de.wikipedia.org/wiki/Aluminium
% atomicWeightAl = 26.9815 * 1.660539040e-27; %u -> kg
% densityAl = 2.7 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t.
% numberConcentrationAl = densityAl / atomicWeightAl;
% molarConcentrationAl = numberConcentrationAl /  jlh.Constants.AvogadroConstant;
% cAl = molarConcentrationAl;

%% concentrations

cAl3p = 0;

% pH = 5.5;
pH = 7; % for electroneutrality in bulk
pOH = 14 - pH;
cH3Op  = 10^(-pH)*10^3; %mol/m^3
cOHm   = 10^(-pOH)*10^3;

cKp = cKNO3_from_kappa(kappa);
cNO3m = cKp; % ignores conductivity contribution of water ions

% cH2O = 55560; % mol/m^3
% cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html
    
% diffusivities
DAl3p = 0.541*10^5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
DH3Op = 9.311*1e-5 / 100^2; % source Haynes 2014, 5-77
DOHm = 5.273*1e-5 / 100^2;  % source Haynes 2014, 5-77
DKp = 1.89*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Arnold 1955 selfdiffusion
% DKp = 1.957*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
% DNO3m = 1700 / 1e6; % mum^2/s = 10^-6 m^2/s, source http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104439&ver=2
DNO3m = 1.902*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77

% charge numbers
zAl3p = -3;
zH3Op = 1;
zOHm = -1;
zKp = 1;
zNO3m = -1;

% environment and physical constants
m.epsilon_r =  jlh.Constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

m.T = 298.15; % K, equivalent to 25 deg C

m.numberOfSpecies = 5; 

m.c_bulk  = [cAl3p,cH3Op,cOHm,cKp,cNO3m]; % bulk concentrations
m.z       = [zAl3p,zH3Op,zOHm,zKp,zNO3m]; % charge number of solute species
m.D       = [DAl3p,DH3Op,DOHm,DKp,DNO3m]; % diffusivity of solute species

m.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
m.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD

% Wbpe  = 0.5e-3/100; % 500mu m
% Wgap  = 0.5e-3/100; %0.5mm / 10000

Wbpe  = 76e-3; %  L= 76 mm
Wgap  = 1e-3; % arbitrary value of 1mm
H  = 0.17e-3/2; % a = 0.17mm, using cell symmetry

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
% E_m = [0.095,0.021,-0.133,-0.286];

% m.PHI_bpe     = E_m(c); % [phi] = V
m.PHI_bulk    = 0;

m.deltaPHI = 1; % arbitrary, 1 V
m.phiSymmetryFactor = 0.5;

%% reactions

% other paramterized reaction parameters
j0a    = 1e-3; % 10^-7 A / cm^2 = 10^-7 * 100^2 A / m^2 = 1e-3 A / m ^2
j0c    = 3e-2; % 3 * 10^-6 A / cm^2 = 3 * 10^-6 * 100^2 A / m^2 = 3e-2 A / m ^2
% ra = na(1-alpha_a) = n(1-beta), rc = nc alpha_c = n beta
% <=> alpha_a = 1 - ra/na
% cases     a   b   c   d   e   f   g
% ra = [1.1e-1,1.1e-1,1.1e-1,1.1e-1,5e-2,5e-2,5e-2]; % fig 3
% rc = [3.5e-2,5.5e-2,7.5e-2,1.1e-1,7.5e-2,1.2e-1,1.7e-2]; % fig 3

ra = 7.5e-2;
rc = 1.1e-1;

m.customParameters.cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
m.customParameters.cH2O = 55560; % mol/m^3

m.customParameters.j0a = j0a;
m.customParameters.j0c = j0c;
m.customParameters.ra = ra;
m.customParameters.rc = rc;

E0a = 1.82; % V, duval2003coupling
E0c = 0.55;
m.customParameters.E0a = E0a;
m.customParameters.E0c = E0c;

%% Acidic Oxygen Reduction, balance to anodic oxygen evolution
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
m.reactions{i}.flux                       = [-1,0,0,0,0]; % no flux

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
m.reactions{i}.flux                       = [0,0,2,0,0]; % inward flux of OHm with cathodic current 

m.nReactions = i;

%% mixed potential
F = jlh.Constants.FaradayConstant;
RT = jlh.Constants.R*jlh.Constants.T;
ja = @(V) j0a*exp(ra*F/RT*(V-E0a));
jc = @(V) -j0c*exp(-rc*F/RT*(V-E0c));
j_tot = @(V) ja(V) + jc(V);

% determine interval
    % mixed potential seems not to always lie between the single standard
    % potentials, especially for irreversible reactions:
%     minE0 = min(E0);
%     maxE0 = max(E0);
minE0 = -2*max(abs([E0a,E0c]));
maxE0 = 2*max(abs([E0a,E0c]));
E_m = fzero( j_tot, [minE0,maxE0]);
m.PHI_bpe = E_m;

%% simulation settings
smootheningFactor = 100;
