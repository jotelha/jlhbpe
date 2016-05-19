% dopeshwarkar2008electrokinetics, small geometry

% concentrations
cTris       = 0.5149; % mM = mol/m^3;
cTrisHp     = 0.4838;
cClm        = 0.4838;
cH3Op       = 7.9e-6;
cOHm        = 1.259e-3;
cBODIPY     = 5e-3;
    
% diffusivities
DTris       = 0.9788e-9; %m^2/s, dhopeshwarkar2008electrokinetics
DTrisHp     = 0.9788e-9; %m^2/s
DClm        = 2.033e-9;
DH3Op       = 5.286e-9;
DOHm        = 5.286e-9;
DBODIPY     = 0.4269e-9;

% charge numbers
zTris = 0;
zTrisHp = 1;
zClm =-1;
zH3Op = 1;
zOHm = -1;
zBODIPY = -2;

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

m.numberOfSpecies = 5; 

m.c_bulk  = [cTrisHp,cClm,cH3Op,cOHm,cBODIPY]; % bulk concentrations
m.z       = [zTrisHp,zClm,zH3Op,zOHm,zBODIPY]; % charge number of solute species
m.D       = [DTrisHp,DClm,DH3Op,DOHm,DBODIPY]; % diffusivity of solute species

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

m.reactions = {};
m.nReactions = 0;