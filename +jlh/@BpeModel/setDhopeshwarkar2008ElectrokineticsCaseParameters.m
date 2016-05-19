function obj = setDhopeshwarkar2008ElectrokineticsCaseParametersSmall(obj)
%% concentrations
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

    obj.bpePoissonEquationBC = 'Robin'; % apply Stern layer Robin BC
    obj.explicitElectrodeGeometry = false;
    %bpePoissonEquationBC = 'Dirichlet'; % apply fixed potential at bpe

    obj.bpeFluxOn = false; % switch species flux due to reactions on or off'

    obj.widthToHeight = 5/3; % width : height relation for single plots
    obj.plotsVisible = 'on'; % 'off' stores plots, but does not show them as pop ups


    %% environment and physical constants
    obj.epsilon_r = constants.RelativePermittivityOfWater; % relative permittivity of electrolyte, here water

    obj.T = 298.15; % K, equivalent to 25 deg C
     
    obj.numberOfSpecies = 5; 

    obj.c_bulk  = [cTrisHp,cClm,cH3Op,cOHm,cBODIPY]; % bulk concentrations
    obj.z       = [zTrisHp,zClm,zH3Op,zOHm,zBODIPY]; % charge number of solute species
    obj.D       = [DTrisHp,DClm,DH3Op,DOHm,DBODIPY]; % diffusivity of solute species

    obj.epsilon     = 0.01; % relation Debye length : simulation domain, lambdaD / L
    obj.delta       = 0.1; % relation width of Stern layer : Debye length, lambdaS / lambdaD
    Wbpe  = 0.5e-3; % 500mu m
    Wgap  = 0.5e-3; %0.5mm / 10000
    
    obj.wMeshRatio = 100; % desired ratio between w_mesh and l (l=1)   
    obj.elementSizeAtSurface = obj.epsilon/10;
    

    obj.Wbpe = Wbpe;
    obj.WinsulatorLeft = 0;
    obj.WinsulatorRight = 0;
    obj.Wcathode = 0;
    obj.Wanode = 0;
    obj.WbulkLeft = Wgap;
    obj.WbulkRight = Wgap;

    obj.PHI_bpe     = 0; % [phi] = V
    obj.PHI_bulk    = 0;
    
    obj.deltaPHI = 7.5;
    obj.phiSymmetryFactor = 0.5;
%     
      %% redox reactions Ox + n e- <=> R
    % k0 = 1e-3; % standard rate constant, 1 cm/s = 1e-3 m/s

    % chemists need to approve concentration boundary conditions
    %% iron ion reduction
%     i = 1;                
%     obj.reactions{i}.name                       = 'aAluminiumOxidation';
%     obj.reactions{i}.reaction                   = 'Al -> Al3+ + 3e-';
%     obj.reactions{i}.oxidants{1}                = 'Al3+';
%     obj.reactions{i}.oxidantConcentrations{1}   = '0'; % no Al3+ in solution
%     obj.reactions{i}.cOx                        = 0; % zero
%     obj.reactions{i}.nuOxidants{1}              = nuAl3p_a;
%     obj.reactions{i}.reductants{1}              = 'Al';
%     obj.reactions{i}.reductantConcentrations{1} = 'cAl';
%     obj.reactions{i}.cR                         = cAl; % absorbed into rate constant
%     obj.reactions{i}.nuReductants{1}            = nuAl_a;
%     obj.reactions{i}.n                          = na;
%     obj.reactions{i}.beta                       = beta_a(c);
%     obj.reactions{i}.k0                         = k0a(c); % mol / (s m^2), after Krawiec 2008
%     obj.reactions{i}.E0                         = E0_a; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
%     obj.reactions{i}.flux                       = [0,0,0,0]; % no flux of H3Op and OHm, should be in the format of -z
% 
%     %%  'Water dissociation'
%     i = i+1;
%     obj.reactions{i}.name                       = 'cWaterReduction'; 
%     % assume only happening in one direction, cathodic
%     obj.reactions{i}.reaction                   = '2H2O + 2e- -> H2 + 2OH-'; % simplified H3O+ + e- <=> H2O + e- -> 1/2 H2 ( + OH- )
% 
%     obj.reactions{i}.oxidants{1}                = 'H2O';
%     % absorbed into rate constant:
%     obj.reactions{i}.oxidantConcentrations{1}   = 'cH2O';
%     obj.reactions{i}.cOx                        = cH2O; % absorbed into rate constant (?)
%     obj.reactions{i}.nuOxidants{1}              = nuH2O_c;
% 
%     obj.reactions{i}.reductants{1}              = 'H2';
%     obj.reactions{i}.reductantConcentrations{1} = '0';
%     obj.reactions{i}.cR                         = 0; % zero
%     obj.reactions{i}.nuReductants{1}            = nuH2_c;
% 
%     % obj.reactions{i}.reductants{2} = 'OH-';
%     % obj.reactions{i}.reductantConcentrations{2} = 'c2';
%     % obj.reactions{i}.nuReductants{2} = 2;
% 
%     obj.reactions{i}.n                          = nc; % as in Krawiec, not sure abouyt it
%     obj.reactions{i}.beta                       = beta_c(c);
%     obj.reactions{i}.k0                         = k0c(c); % m^4 / (s mol)
%     obj.reactions{i}.E0                         = E0_c; % V, source: https://en.wikipedia.org/wiki/Standard_electrode_potential_(data_page)
%     %obj.reactions{i}.flux = [0,1]; % inward flux of OHm with cathodic current 
%     obj.reactions{i}.flux                       = [0,0,0,0]; % inward flux of OHm with cathodic current 
% 
%     obj.nReactions = i;
    obj.reactions = {};
    obj.nReactions = 0;
end