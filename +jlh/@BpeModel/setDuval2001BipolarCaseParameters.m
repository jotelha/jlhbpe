function obj = setDuval2001BipolarCaseParameters(obj,c)
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

%% anodic reaction
na          = 3;
nuAl_a      = 1; % reductant
nuAl3p_a    = -1; % oxidant
E0_a        = -1.82; %V, after Duval in pH = 5.5
% E0_a = -1.662; %V, after https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page), ref 1

%% cathodic reaction
nc          = 2;
nuOHm_c     = 2;
nuH2_c      = 1;
nuH2O_c     = -2;
E0_c        = -0.55; %V, after Duval
% E0_c      = -0.8277; %V, after https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page), ref 7
%% other paramterized reaction parameters
i0a    = [0.08,0.10,0.13,5.32] / 100; % muA /cm^2 = 100^2 / m^2 * A / 1000^2 = 0.01 A /m^2
i0c    = [2.12,1.57,1.16,0.21] / 100;
% ra = na(1-alpha_a) = n(1-beta), rc = nc alpha_c = n beta
% <=> alpha_a = 1 - ra/na
ra = [7,6.5,5.6,1.6]/100;
rc = [10.7,10.3,9.3,6.9]/100;
beta_a = 1 - ra/na;
beta_c = rc/nc;


E_m = [-0.58,-0.62,-0.65,-1.76];

%% concentrations
cKNO3  = [100,50,10,0.2]; % mmol /l = mol/m^3

% pure aluminium: https://de.wikipedia.org/wiki/Aluminium
atomicWeightAl = 26.9815 * 1.660539040e-27; %u -> kg
densityAl = 2.7 * 100^3/1000; %g/cm^3 -> kg/m^3 at r.t.
numberConcentrationAl = densityAl / atomicWeightAl;
molarConcentrationAl = numberConcentrationAl / constants.AvogadroConstant;
cAl = molarConcentrationAl;

pH = 5.5;
pOH = 14 - pH;
cH3Op  = 10^(-pH)*10^3; %mol/m^3
cOHm   = 10^(-pOH)*10^3;

cKp = cKNO3(c);
cNO3m = cKNO3(c);

cH2O = 55560; % mol/m^3
cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html

%% kinetics
% i0 = n F k0 cOx^(1-beta) cR^beta
% k0 = i0/(n F*cOx^(1-beta) * cR^beta) = i0/(n F c_ref)
% irreversible reactions
cref_a = cAl;
cref_c = cH2O;
k0a = i0a ./ (na*jlh.Constants.FaradayConstant*cref_a);
k0c = i0c ./ (nc*jlh.Constants.FaradayConstant*cref_c);

% 2H2O + 2e- -> H2 + 2OH-
% after Newman: sum si Mi = n e-
% sH2 = 1; sOHm = 2; sH2O = -2;
% i0 = n F k0 prod ci^(qi+beta*si)
%   si > 0: pi = si, qi = 0
%   si < 0: pi = 0, qi = -si
% i0 = n F k0 cH2^beta cOHm^(2*beta) cH2O^(2*(1-beta))
% cOx = cH2O^2
% cR = cH2 * cOHm^2
k0c_rev = i0c ./ (jlh.Constants.FaradayConstant .* ...
        cH2.^(abs(nuH2_c).*beta_c) .* ...
        cOHm.^(abs(nuOHm_c).*beta_c) .* ...
        cH2O.^(abs(nuH2O_c).*(1-beta_c)));   
    
%% diffusivities
DKp = 1.89*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Arnold 1955 selfdiffusion
% DKp = 1.957*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
% DNO3m = 1700 / 1e6; % mum^2/s = 10^-6 m^2/s, source http://bionumbers.hms.harvard.edu/bionumber.aspx?&id=104439&ver=2
DNO3m = 1.902*1e-5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
% DAl3p = 0.541*10^5 / 100^2; % cm^2/s = 100^-2 m^2/s, source Haynes 2014, 5-77
DH3Op = 9.311*1e-5 / 100^2; % source Haynes 2014, 5-77
DOHm = 5.273*1e-5 / 100^2;  % source Haynes 2014, 5-77

%% electrical conductivity:
% KNO3 in aqueous solution
xiKNO3_haynes   = [0.5, 1, 2, 5, 10, 15, 20, 25];
% kappa in mS/cm concentration in mass percent
%   c       0.5     1       2       5       10      15      20      25
%   kappa   5.5     10.7    20.1    47.0    87.3    124     157     182
% source: Electrical conductivity of aqueous solutions
%   Hayne 2014 CRC Handbook of chemistry and physics, 5-73

cKNO3_haynes    = [0.05, 0.099, 0.2, 0.509, 1.051, (1.509+1.747)/2, 2.24, 2.759] * 10^3; % mol / L = 10^3 mol/m^3
% source: Electrical conductivity of aqueous solutions
%   Hayne 2014 CRC Handbook of chemistry and physics, 5-73
kappa_haynes = [5.5,10.7,20.1,47,87.3,124,157,182] *0.1; %mS/cm = 100/1000 S/m

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
coeffs = polyfit(cKNO3_haynes_dilute,kappa_haynes_dilute,2);
kappa = polyval(coeffs,cKNO3);

%% switches

% Duval formulation:    i = i0 * (exp( ra f (E-E0) ) - exp( - rc f (E-E0) )
% Newman form:          i = i0 * (exp( (1-beta) n f eta ) - exp( -beta n f eta)
% Bard form             i = n F k0 ( cOx*exp(-beta f 

    obj.bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
    obj.explicitElectrodeGeometry = false;
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
    obj.customParameters.cAl = cAl; 
    obj.customParameters.cH2O = cH2O; % mol/m^3
%     obj.customParameters.cO2 = 258e-3; % 258 muMol = 258e-6 mol/l = 258e-3 mol/m^3? http://www.engineeringtoolbox.com/oxygen-solubility-water-d_841.html
    %cO2 = 0.29  % 0.29mM = 0.29e-3 mol/l = 0.29 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html
    % cH2 = 0; % assume no hydrogen in solution
%     obj.customParameters.cH2 = 0.44e-6; % 0.44nM = 0.44e-9 M = 0.44e-6 mol/m^3, http://www1.lsbu.ac.uk/water/electrolysis.html

    obj.numberOfSpecies = 4; %H3O+,OH-,K+,NO3-
%     pH = 7;
%     cH3Op_bulk  = 10^(-pH)*10^3; %mol/m^3
%     cOHm_bulk   = cH3Op_bulk;
%     cFe2p_bulk  = 0;

    %obj.customParameters.cH3Op   = cH3Op;
    %obj.customParameters.cOHm    = cOHm;
    %obj.customParameters.cKp    = cOHm;
    obj.customParameters.cAl   = cAl;
    obj.customParameters.cH2o = cH2O;


    % H3O+,OH-,Na+,DS-,ClO4-
%     obj.c_bulk  = [cH3Op_bulk,cOHm_bulk,cNap_bulk,cDSm_bulk,cClO4m_bulk]; % bulk concentrations
%     obj.z       = [1,-1,1,-1,-1]; % charge number of solute species
%     obj.D       = [jlh.Constants.DiffusivityH3Op,1e-9,1e-9,1e-9,1e-9]; % diffusivity of solute species
    obj.c_bulk  = [cH3Op,cOHm,cKp,cNO3m]; % bulk concentrations
    obj.z       = [1,-1,1,-1]; % charge number of solute species
    obj.D       = [DH3Op,DOHm,DKp,DNO3m]; % diffusivity of solute species
% 
  
    
    obj.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
    obj.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
%     obj.gamma       = obj.lambdaD/W; % relation Debye length : simulation domain width, lambdaD / W0
%     obj.beta        = obj.lambdaD/Wbpe;

    L     = obj.lambdaD/obj.epsilon;

    % measures 
    a = [0.15,0.20,2.3] / 1000; % mm / 1000 = m, the gap between two plates
    L0 = 76 / 1000; %mm / 1000 = m, the length of BPE
    
    
    % Wbpe = L0;
    Wbpe  = 8e-5; % 10mm / 10000
    Wgap  = 1/10*Wbpe; %0.5mm / 10000
    
    obj.wMeshRatio = 100; % desired ratio between w_mesh and l (l=1)   
    obj.elementSizeAtSurface = obj.epsilon/10;
    

%     obj.Wbpe = Wbpe;
%     obj.WinsulatorLeft = Wgap;
%     obj.WinsulatorRight = Wgap;
%     obj.Wcathode = Wel;
%     obj.Wanode = Wel;
%     obj.WbulkLeft = L;
%     obj.WbulkRight = L;
    obj.Wbpe = Wbpe;
    obj.WinsulatorLeft = 0;
    obj.WinsulatorRight = 0;
    obj.Wcathode = 0;
    obj.Wanode = 0;
    obj.WbulkLeft = Wgap;
    obj.WbulkRight = Wgap;

    obj.PHI_bpe     = E_m(c); % [phi] = V
    obj.PHI_bulk    = 0;
    
    obj.deltaPHI = 1*(1+2*Wgap/Wbpe);
    obj.phiSymmetryFactor = 0.5;
%     
      %% redox reactions Ox + n e- <=> R
    % k0 = 1e-3; % standard rate constant, 1 cm/s = 1e-3 m/s

    % chemists need to approve concentration boundary conditions
    %% iron ion reduction
    i = 1;                
    obj.reactions{i}.name                       = 'aAluminiumOxidation';
    obj.reactions{i}.reaction                   = 'Al -> Al3+ + 3e-';
    obj.reactions{i}.oxidants{1}                = 'Al3+';
    obj.reactions{i}.oxidantConcentrations{1}   = '0'; % no Al3+ in solution
    obj.reactions{i}.cOx                        = 0; % zero
    obj.reactions{i}.nuOxidants{1}              = nuAl3p_a;
    obj.reactions{i}.reductants{1}              = 'Al';
    obj.reactions{i}.reductantConcentrations{1} = 'cAl';
    obj.reactions{i}.cR                         = cAl; % absorbed into rate constant
    obj.reactions{i}.nuReductants{1}            = nuAl_a;
    obj.reactions{i}.n                          = na;
    obj.reactions{i}.beta                       = beta_a(c);
    obj.reactions{i}.k0                         = k0a(c); % mol / (s m^2), after Krawiec 2008
    obj.reactions{i}.E0                         = E0_a; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    obj.reactions{i}.flux                       = [0,0,0,0]; % no flux of H3Op and OHm, should be in the format of -z

    %%  'Water dissociation'
    i = i+1;
    obj.reactions{i}.name                       = 'cWaterReduction'; 
    % assume only happening in one direction, cathodic
    obj.reactions{i}.reaction                   = '2H2O + 2e- -> H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )

    obj.reactions{i}.oxidants{1}                = 'H2O';
    % absorbed into rate constant:
    obj.reactions{i}.oxidantConcentrations{1}   = 'cH2O';
    obj.reactions{i}.cOx                        = cH2O; % absorbed into rate constant (?)
    obj.reactions{i}.nuOxidants{1}              = nuH2O_c;

    obj.reactions{i}.reductants{1}              = 'H2';
    obj.reactions{i}.reductantConcentrations{1} = '0';
    obj.reactions{i}.cR                         = 0; % zero
    obj.reactions{i}.nuReductants{1}            = nuH2_c;

    % obj.reactions{i}.reductants{2} = 'OH-';
    % obj.reactions{i}.reductantConcentrations{2} = 'c2';
    % obj.reactions{i}.nuReductants{2} = 2;

    obj.reactions{i}.n                          = nc; % as in Krawiec, not sure abouyt it
    obj.reactions{i}.beta                       = beta_c(c);
    obj.reactions{i}.k0                         = k0c(c); % m^4 / (s mol)
    obj.reactions{i}.E0                         = E0_c; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
    %obj.reactions{i}.flux = [0,1]; % inward flux of OHm with cathodic current 
    obj.reactions{i}.flux                       = [0,0,0,0]; % inward flux of OHm with cathodic current 

    obj.nReactions = i;
end