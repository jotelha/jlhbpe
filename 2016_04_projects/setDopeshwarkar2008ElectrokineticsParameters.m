% dopeshwarkar2008electrokinetics, small geometry

% concentrations
cTris       = 0.5149; % mM = mol/m^3;
cTrisHp     = 0.4838;
cClm        = 0.4838;
cH3Op       = 7.9e-6;
cOHm        = 1.259e-3;
cBODIPY     = 5e-3;

%% diffusivities
DTris       = 0.9788e-9; %m^2/s, dhopeshwarkar2008electrokinetics
DTrisHp     = 0.9788e-9; %m^2/s
DClm        = 2.033e-9;
DH3Op       = 5.286e-9;
DOHm        = 5.286e-9;
DBODIPY     = 0.4269e-9;

%% charge numbers
zTris = 0;
zTrisHp = 1;
zClm =-1;
zH3Op = 1;
zOHm = -1;
zBODIPY = -2;

%% switches

% Duval formulation:    i = i0 * (exp( ra f (E-E0) ) - exp( - rc f (E-E0) )
% Newman form:          i = i0 * (exp( (1-beta) n f eta ) - exp( -beta n f eta)
% Bard form             i = n F k0 ( cOx*exp(-beta f 

m.bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
m.explicitElectrodeGeometry = false;
%bpePoissonEquationBC = 'Dirichlet'; % apply fixed potential at bpe

m.bpeFluxOn = false; % switch species flux due to reactions on or off'

m.widthToHeight = 5/3; % width : height relation for single plots
m.plotsVisible = 'on'; % 'off' stores plots, but does not show them as pop ups


%% environment and physical constants
m.epsilon_r = constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

m.T = 298.15; % K, equivalent to 25 deg C

m.numberOfSpecies = 5; 

m.c_bulk  = [cTrisHp,cClm,cH3Op,cOHm,cBODIPY]; % bulk concentrations
m.z       = [zTrisHp,zClm,zH3Op,zOHm,zBODIPY]; % charge number of solute species
m.D       = [DTrisHp,DClm,DH3Op,DOHm,DBODIPY]; % diffusivity of solute species

m.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
m.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
Wbpe  = 0.5e-3; % 500mu m
Wgap  = 0.5e-3; %0.5mm / 10000

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

m.deltaPHI = 7.5;
m.phiSymmetryFactor = 0.5;
%     
  %% redox reactions Ox + n e- <=> R
% k0 = 1e-3; % standard rate constant, 1 cm/s = 1e-3 m/s

% chemists need to approve concentration boundary conditions
%% iron ion reduction
%     i = 1;                
%     m.reactions{i}.name                       = 'aAluminiumOxidation';
%     m.reactions{i}.reaction                   = 'Al -> Al3+ + 3e-';
%     m.reactions{i}.oxidants{1}                = 'Al3+';
%     m.reactions{i}.oxidantConcentrations{1}   = '0'; % no Al3+ in solution
%     m.reactions{i}.cOx                        = 0; % zero
%     m.reactions{i}.nuOxidants{1}              = nuAl3p_a;
%     m.reactions{i}.reductants{1}              = 'Al';
%     m.reactions{i}.reductantConcentrations{1} = 'cAl';
%     m.reactions{i}.cR                         = cAl; % absorbed into rate constant
%     m.reactions{i}.nuReductants{1}            = nuAl_a;
%     m.reactions{i}.n                          = na;
%     m.reactions{i}.beta                       = beta_a(c);
%     m.reactions{i}.k0                         = k0a(c); % mol / (s m^2), after Krawiec 2008
%     m.reactions{i}.E0                         = E0_a; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
%     m.reactions{i}.flux                       = [0,0,0,0]; % no flux of H3Op and OHm, should be in the format of -z
% 
%     %%  'Water dissociation'
%     i = i+1;
%     m.reactions{i}.name                       = 'cWaterReduction'; 
%     % assume only happening in one direction, cathodic
%     m.reactions{i}.reaction                   = '2H2O + 2e- -> H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )
% 
%     m.reactions{i}.oxidants{1}                = 'H2O';
%     % absorbed into rate constant:
%     m.reactions{i}.oxidantConcentrations{1}   = 'cH2O';
%     m.reactions{i}.cOx                        = cH2O; % absorbed into rate constant (?)
%     m.reactions{i}.nuOxidants{1}              = nuH2O_c;
% 
%     m.reactions{i}.reductants{1}              = 'H2';
%     m.reactions{i}.reductantConcentrations{1} = '0';
%     m.reactions{i}.cR                         = 0; % zero
%     m.reactions{i}.nuReductants{1}            = nuH2_c;
% 
%     % m.reactions{i}.reductants{2} = 'OH-';
%     % m.reactions{i}.reductantConcentrations{2} = 'c2';
%     % m.reactions{i}.nuReductants{2} = 2;
% 
%     m.reactions{i}.n                          = nc; % as in Krawiec, not sure abouyt it
%     m.reactions{i}.beta                       = beta_c(c);
%     m.reactions{i}.k0                         = k0c(c); % m^4 / (s mol)
%     m.reactions{i}.E0                         = E0_c; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
%     %m.reactions{i}.flux = [0,1]; % inward flux of OHm with cathodic current 
%     m.reactions{i}.flux                       = [0,0,0,0]; % inward flux of OHm with cathodic current 
% 
%     m.nReactions = i;
m.reactions = {};
m.nReactions = 0;