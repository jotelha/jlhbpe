% 
% all model settings except create commands
%
function obj = updateModel(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    fprintf('  Modify settings of model "%s"...\n',obj.model_tag);

    % ModelUtil.showProgress(true);               % shows computation progress in MATLAB

    %% note concerning comsol normal vectors
    % source: https://www.comsol.de/community/forums/general/thread/9119/
    % the normal direction are defined wrt the body/domain, it points normally
    % "outward" from a closed volume ("up") and if you define a plane you need
    % to say where is the geoemtrical entity of higher value (volume for a
    % plane, surface for an edge, but a point has a singularity w.r.t normal ;)
    % to define the outward ("up") direction ofthe normal.


    %% prepare replacement patterns
    % %u or %d for integers, %s for strings
    % all identifiers have to be single column arrays

    % pdeConvectionCoefficientAlpha = 'z%d*D%d/RT*F_const*phix';
    % pdeConvectionCoefficientAlpha = sprintf('z%1$d*D%1$d/RT*F_const*phix',1:numberOfSpecies);
    % convection coefficient alpha = zi * Di * R /RT * phix

    % derived expressions
    lengthL     = 'lambdaD/epsilon'; % domain depth
    widthW      = 'lambdaD/gamma'; % domain width
    widthWbpe   = 'lambdaD/beta'; 
    lambdaS     = 'delta*lambdaD';
    ionicStrengthSummand = prepTerm('z_id^2*c_bulk_id','z_id','c_bulk_id',obj.z_id,obj.c_bulk_id);
    ionicStrengthTerm = strcat('1/2*(', strjoin( ionicStrengthSummand','+' ), ')');
    debyeLengthTerm = 'sqrt(epsilon0_const*epsilon_r*RT/(2*F_const^2*IonicStrength))';

    % related to reactions
    %fluxTerm = cell(nReactions,1);
    % reaction.flux actiually is each species stochiometric coefficient nu for
    % this reaction Ox + n e -> R
    %   nu > 0: oxidant
    %   nu < 0: reductant
    fluxDirections = cellfun( @(r) r.flux', obj.reactions, 'UniformOutput',false);
    fluxDirections = [fluxDirections{:}];
    fluxSummandPattern = prepTerm('-(nu_id*i_id/(F_const*n_id))','i_id','n_id',obj.i_id,obj.n_id);
    specificFluxSummand = arrayfun( @(n) prepTerm(fluxSummandPattern{n},'nu_id',obj.nu_id(:,n)),1:obj.nReactions,'UniformOutput',false);
    specificFluxSummand = [specificFluxSummand{:}];

    fluxTerm = cell(obj.numberOfSpecies,1);
    % positiveFlux = cell(numberOfSpecies,1);
    % negativeFlux = cell(numberOfSpecies,1);
    for i=1:obj.numberOfSpecies 
    %     positiveFluxSummand = specificFluxSummand{i,fluxDirections(i,:) > 0};
    %     negativeFluxSummand = specificFluxSummand{i,fluxDirections(i,:) < 0};
    %     positiveFlux{i} = strjoin(positiveFluxSummand','+');
    %     negativeFlux{i} = strjoin( strcat('-',negativeFluxSummand)','');
    %     fluxTerm{i} = strcat(positiveFlux{i},negativeFlux{i});
        nonZeroFluxSummand = specificFluxSummand(i,fluxDirections(i,:)~=0);
        fluxTerm{i} = strjoin( nonZeroFluxSummand, '' );
        if strcmp(fluxTerm{i},'')
            fluxTerm{i} = '0';
        end
    end
    % fluxTerm{ strcmp(fluxTerm,'') } = '0';

    % concentrations of species participating in surface reactions
    oxidantConcentrations = cellfun( @(r) r.oxidantConcentrations{1}, obj.reactions','UniformOutput',false);
    reductantConcentrations = cellfun( @(r) r.reductantConcentrations{1}, obj.reactions','UniformOutput',false);

    % IMPORTANT 
    % opposite of convention in comsol
    %   anodic currents NEGATIVE, cathodic currents POSITIVE

    % comsol convention:
    %   anodic currents POSITIVE, cathodic currents NEGATIVE

    % anodic currents    
    anodicTerms = prepTerm(...
        'n_id*F_const*k0_id*(c_R*exp((1-beta_id)*F_const/RT*(V-E0_id)))',...
        'n_id','k0_id','c_R','beta_id','E0_id',...
        obj.n_id,obj.k0_id,reductantConcentrations,obj.beta_id,obj.E0_id);
    anodicCurrent = strjoin(anodicTerms','+');

    % cathodic currents
    cathodicTerms = prepTerm(...
        '-n_id*F_const*k0_id*(c_Ox*exp(-beta_id*F_const/RT*(V-E0_id)))',...
        'n_id','k0_id','c_Ox','beta_id','E0_id',...
        obj.n_id,obj.k0_id,oxidantConcentrations,obj.beta_id,obj.E0_id);    
    cathodicCurrent = strjoin(cathodicTerms','');

    totalCurrent = strjoin({anodicCurrent,cathodicCurrent},'');

    % current per reaction
    % Butler-Volmer, cathodic direction > 0 in form of, as in Bard p. 96
    % i = n*F/RT*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
    % in anodic direction > 0:
    reactionCurrent = prepTerm(...
        'n_id*F_const*k0_id*( - c_Ox*exp(-beta_id*F_const/RT*(V-E0_id)) + c_R*exp((1-beta_id)*F_const/RT*(V-E0_id)) ) ',...
        'n_id','k0_id','c_Ox','c_R','beta_id','E0_id',...
        obj.n_id,obj.k0_id,oxidantConcentrations,reductantConcentrations,obj.beta_id,obj.E0_id);
    %         n_id{i},'*F_const*',k0_id{i},'*( ',...
    %         reactions{i}.oxidantConcentrations{1},...
    %         '*exp(-',beta_id{i},'*F_const/RT * (V-',E0_id{i},'))',...
    %         '-',...
    %         reactions{i}.reductantConcentrations{1},...
    %         '*exp((1-',beta_id{i},')*F_const/RT * (V-',E0_id{i},')) )'...
    %         ), ['current due to surface reaction ', reactions{i}.reaction]);


    %% original formulation 
    % in comsol, coefficient steady state PDE have the form
    %   div*(-c*grad u - alpha*u + gamma) + beta*grad u = f 
    % we solve
    %   div*(-D*grad c - z*D*F/RT*phix*c ) = 0 (NP - Nernst-Planck)
    %   div*(-epsilon_r*epsilon_0*grad phi) = F*sum(z_i*c_i) (P - Poisson)
    % thus for NP
    %   diffusion coefficient c = D
    %   conservative flux convection coefficient alpha = z*D*F/RT*phix
    %   source term f = 0
    % and for P
    %   diffusion coefficient c = epsilon_r*epsilon_0
    %   source term f = F*sum(z_i*c_i)

    %% dimensionless formulation
    % we solve
    %   div*(-grad c - z*phix*c ) = 0 (NP - Nernst-Planck)
    %   div*(-epsilon^2*grad phi) = 1/2*sum(z_i*c_i) (P - Poisson)
    % thus for NP
    %   diffusion coefficient c = 1
    %   conservative flux convection coefficient alpha = z*phix
    %   source term f = 0
    % and for P
    %   diffusion coefficient c = epsilon^2
    %   source term f = 1/2*sum(z_i*c_i)
    % here, epsilon = lambdaD/L with
    %   Debye length lambdaD = sqrt(epsilon_r*epsilon_0*RT/2F^2I)
    %   ionic strength I = 1/2 * sum( z_i^2 c_i_bulk)

    % pdeConvectionCoefficientAlpha = prepTerm('z_id*D_id*F_const/RT*phix','z_id','D_id',z_id,D_id);

    %% Nernst-Planck equation 
    % domain terms NP:
    pdeNPDiffusionCoefficientC = '1';
    pdeNPConvectionCoefficientAlphaX = prepTerm('z_id*phix','z_id',obj.z_id);
    pdeNPConvectionCoefficientAlphaY = prepTerm('z_id*phiy','z_id',obj.z_id);

    pdeNPSourceTerm = strzeros(obj.numberOfSpecies,1); % dropped '

    % for log concentration PDE
    % pdeGammaTerm = prepTerm('-(D_id*C_idx+z_id*D_id*F_const/RT*phix)',...
    %     'D_id','C_id','z_id',D_id,C_id,z_id); % minus convention
    % pdeSourceTerm = prepTerm('C_idx*(D_id*C_idx+z_id*D_id*F_const/RT*phix)',...
    %     'D_id','C_id','z_id',D_id,C_id,z_id);

    % boundary terms NP:
    pdeNPbulkDirichletBC    = prepTerm('c_bulk_id/c_ref','c_bulk_id',obj.c_bulk_id); % fixed concentration at bulk
    pdeNPbpeNeumannBC       = prepTerm('L*N_id/(D_id*c_ref)','N_id','D_id',obj.N_id,obj.D_id); % species flux at surface

    c_ref = 'IonicStrength';

    % initial value terms
    pdeNPinitialValues = prepTerm('c_pb_id(y,phi_bpe)','c_pb_id',obj.c_pb_id);
    %% Poisson equation
    pdePDiffusionCoefficientC = 'epsilon^2';
    % chargeDensitySummand = prepTerm('z_id*c_id','z_id','c_id',z_id,c_id);
    % chargeDensityTerm = strcat('F_const*(',...
    %     strjoin(chargeDensitySummand','+' ), ')');
    chargeDensitySummand = prepTerm('z_id*c_id','z_id','c_id',obj.z_id,obj.c_id);
    chargeDensityTerm = strcat('1/2*(',... % 1/2 arises in the chosen dimensionless formulation
        strjoin(chargeDensitySummand','+' ), ')');

    pdePSourceTerm = prepTerm('chargeDensityRampFactor*chargeDensityTerm','chargeDensityTerm',chargeDensityTerm); % allow for ramping


    % boundary terms P

    RT_term = 'R_const*T'; 
    UT_term = 'RT/F_const'; % thermal voltage
    %phi0_term = 'PHI_bpe/UT'; % dimensionless surface potential

    % pdePbulkDirichletBC   = 'phi_bulk/UT'; % voltage at bulk, usually zero
%     pdePbulkDirichletBC   = 'phi_we-x*L/W*deltaPhi'; % voltage at bulk, usually zero
    pdePbulkDirichletBC   = 'phi_pb(1,phi_we)-x*L/W*(phi_pb(1,phi_we)-phi_pb(1,phi_ce))'; % voltage at bulk, usually zero

    wePotentialTerm       = 'phi_bpe+phiSymmetryFactor*deltaPhi';
    cePotentialTerm       = 'phi_bpe-(1-phiSymmetryFactor)*deltaPhi';
    
    deltaPhiTerm          = sprintf('deltaPhiRampFactor*%e',obj.delta_phi);


    % robin bc in comsol are of form
    %   -n*(-c grad u - alpha u + gamma) = g - q*u
    % with dimensionless phi0 = PHI0 / UT we have (UT = RT/F):
    %   phi + lambdaS/L*n*grad phi = phi0, or rearranged
    %   -n*(-epsilon^2 grad phi) = epsilon^2*L/lambdaS*phi0 - epsilon^2*L/lambdaS*phi
    % with capacitcance of Stern layer C_Stern = epsilon0_const*epsilon_r/lambdaS
    % and its dimensionless equivalent 
    %   C*_Stern = epsilon^2*L/lambdaS = epsilon^2 *(L/lambdaD)*(lambdaD/lambdaS) 
    %            = epsilon/delta, we have
    %   -n*(-epsilon^2 grad phi) = C*_Stern*phi0 - C*_Stern*phi
    % notice:
    %   epsilon^2 = lambdaD^2/L^2 => C*_Stern = lambdaD^2/(lambdaS*L)
    %   CStern = 'epsilon0_const*epsilon_r/lambdaS';
    % we identify
    %   boundary flux g = C*_Stern*phi0/UT
    %   boundary adsorption term q = C*_Stern
    % C_Stern = 'epsilon^2*L/lambdaS'; % dimensionless Stern layer capacitance
    C_Stern = 'epsilon/delta'; 
    C_Stern_dimensional = 'epsilon0_const*epsilon_r/(epsilon^2*L)*C_Stern';

    pdePbpeRobinBCg     = 'C_Stern*phi_s';
    pdePbpeRobinBCq     = 'C_Stern';

    % for delta = 0 => lambdaS = 0, BC reduces to
    %   phi = phi0

    pdePbpeDirichletBC = 'phi_s';

    % for working and counter electrodes
    pdePelectrodeDirichletBC = 'phi_pb(y,phi_s)';

    % initial values
    pdePinitialValues = 'phi_pb(y,phi_bpe)';

    %% analytic Poisson-Boltzmann solution, uses concentration of first species
    % assuming symmetric z:z electrolyte
    % phiPBFunction = '4*RT/(abs(z1)*F_const)*atanh( tanh( abs(z1)*F_const * phi_bpe/(4*RT) )*exp(-x/lambdaD))';
    % phiPBxFunction = '- 4*RT/(lambdaD*abs(z1)*F_const)* sinh( abs(z1) * F_const* phi_bpe/(4*RT) )* cosh( abs(z1) * F_const* phi_bpe/(4*RT) ) / ( exp(x/lambdaD)* cosh( abs(z1) * F_const* phi_bpe/(4*RT) )^2 - exp(-x/lambdaD) * sinh( abs(z1) * F_const* phi_bpe/(4*RT) )^2)';
    % correction for Stern layer:
    % phiPBFunction = '4*RT/(abs(z1)*F_const)*atanh( tanh( abs(z1)*F_const * phi_bpe/(4*RT) )*exp(-(x+lambdaS)/lambdaD))';
    % dimensionless: 
    %   phiPB* = phiPB/UT = F/RT*phiPB
    %   phiPBx* = L/UT * phiPBx
    phiPBFunction = '4/abs(z1)*atanh( tanh( abs(z1)*phi0/4 )*exp(-(x/epsilon+delta)))'; % dimensionless
    % phiPBxFunction = '- 4*RT/(lambdaD*abs(z1)*F_const)* sinh( abs(z1) * F_const* phi_bpe/(4*RT) )* cosh( abs(z1) * F_const* phi_bpe/(4*RT) ) / ( exp((x+lambdaS)/lambdaD)* cosh( abs(z1) * F_const* phi_bpe/(4*RT) )^2 - exp(-(x+lambdaS)/lambdaD) * sinh( abs(z1) * F_const* phi_bpe/(4*RT) )^2)';
    phiPBxFunction = '- 4/(epsilon*abs(z1))* sinh( abs(z1) * phi0/4 )* cosh( abs(z1) * phi0/4 ) / ( exp(x/epsilon+delta)* cosh( abs(z1) * phi0/4 )^2 - exp(-(x/epsilon+delta)) * sinh( abs(z1) * phi0/4 )^2)';


    % cPBFunction = prepTerm('c_bulk_id*exp(-z_id*F_const*phi_pb(x)/RT)',...
    %     'c_bulk_id','z_id',c_bulk_id,z_id);
    % dimensionless:
    %   cPB* = c_bulk/c_ref * exp( -z*phiPB*)
    cPBFunction = prepTerm('c_bulk_id/c_ref*exp(-z_id*phi_pb(x,phi0))',...
        'c_bulk_id','z_id',obj.c_bulk_id,obj.z_id);

    % for i=1:nReactions 
    %     %fluxDirections = reactions{i}.flux';
    %     % each reaction might have several species reacting
    %     % a single species flux looks like N = +/- i / (n*F)
    %     % here i is the current density due to surface reaction
    %     % n is the number of electrons participating (?)
    %     % the sign depends on whether species is participating as reductant or
    %     % oxidant
    %     fluxSummandPattern = prepTerm('(flux_id*i_id/(F_const*n_id))','i_id','n_id',i_id{i},n_id{i});
    %     fluxSummand{i} = prepTerm(fluxSummandPattern{1},'flux_id',cellstr(int2str(fluxDirections)));    
    %     %activeFluxSummand{i} = (fluxDirections ~= 0);
    %     %fluxTerm{i} = strcat(strjoin( activeFluxSummand', '+'));
    % %     if isempty(fluxTerm{i})
    % %         fluxTerm{i} = '0'; % zero flux 
    % %     else
    % %         fluxTerm{i} = strcat('(',fluxTerm{i},')');
    % %     end
    % end  

    %% set gloabl parameters
    fprintf('  Setting global parameters...\n');

    % for general pde
    % gammaTerm = '(D%i*C%ix+z%i*D%i*F_const/RT*phix)';
    % sourceTerm = 'C%ix*(D%i*C%ix+z%i*D%i*F_const/RT*phix)';

    for n=1:obj.numberOfSpecies
        obj.m.param.set(obj.D_id{n}, attachUnit(obj.D(n),'m^2/s'), 'Diffusion coefficient');
        obj.m.param.set(obj.z_id{n}, attachUnit(obj.z(n),'1'), 'Charge number');
        obj.m.param.set(obj.c_bulk_id{n}, attachUnit(obj.c_bulk(n),'mol/m^3'), 'Bulk concentration');
    end

    % all below should be changed to directly insert numeric values from matlab
    obj.m.param.set('surfacePotentialRampFactor', 1);
    obj.m.param.set('surfaceFluxRampFactor', 1);
    obj.m.param.set('chargeDensityRampFactor', 1);
    obj.m.param.set('deltaPhiRampFactor', 1);

    
    obj.m.param.set('epsilon', num2str(obj.epsilon), 'Dimensionless Debye length scale');
    obj.m.param.set('delta', num2str(obj.delta), 'Dimensionless Stern layer width');
    obj.m.param.set('gamma', num2str(obj.gamma), 'Dimensionless relation betweeb Debye length and domain width');
    obj.m.param.set('beta', num2str(obj.beta), 'Dimensionless relation betweeb Debye length and BPE width');

    obj.m.param.set('W', widthW, 'Cell width');
    obj.m.param.set('Wbpe', widthWbpe, 'Width of BPE');
    obj.m.param.set('L', lengthL, 'Cell depth');
    
    obj.m.param.set('lambdaD', debyeLengthTerm, 'Debye length');
    obj.m.param.set('kappa', 'lambdaD^(-1)', 'reciprocal of Debye length');
    obj.m.param.set('lambdaS', lambdaS); % width of Stern layer, move to parameter settings

    obj.m.param.set('T', attachUnit(obj.T,'K'), 'Temperature'); % move to parameter settings
    obj.m.param.set('RT', RT_term, 'Molar gas constant * Temperature');
    obj.m.param.set('UT', UT_term, 'thermal voltage, R*T/F');
    obj.m.param.set('epsilon_r', num2str(obj.epsilon_r), 'relative permittivity of electrolyte');

    obj.m.param.set('IonicStrength',ionicStrengthTerm, '');
    obj.m.param.set('c_ref',c_ref, 'reference concentration for dimensionless formulation');

    % m.param.set('C_Stern', 'epsilon0_const*epsilon_r/lambdaS', 'Stern layer capacitance'); % Stern layer capacitance
    obj.m.param.set('C_Stern', C_Stern, 'Stern layer capacitance, dimensionless'); % Stern layer capacitance
    obj.m.param.set('C_Stern_dimensional', C_Stern_dimensional, 'Stern layer capacitance'); % Stern layer capacitance

    % m.param.set('phi_bpe',attachUnit(phi_bpe,'V'),'potential at surface of bpe');
    obj.m.param.set('PHI_bpe',['surfacePotentialRampFactor*',attachUnit(obj.PHI_bpe,'V')],'potential at surface of bpe');
    obj.m.param.set('phi_bpe',obj.phi_bpe,'dimensionless potential at surface of bpe');
    obj.m.param.set('deltaPhi',deltaPhiTerm,'potential difference between working electrode and counter electrode');
    obj.m.param.set('phiSymmetryFactor',obj.phiSymmetryFactor,'symmetry of applied potential around bpe potential');

    obj.m.param.set('phi_we',wePotentialTerm,'dimensionless potential at working electrode');
    obj.m.param.set('phi_ce',cePotentialTerm,'dimensionless potential at counter electrode');

    obj.m.param.set('phi_bulk','0','bulk electrolyte as reference');

    % concentration of water, iron and iron ions
    obj.m.param.set('cH2O', attachUnit(obj.cH2O,'mol/m^3'), '1000 g / l * 1 mol / 18 g = 55.56 mol / L');
    obj.m.param.set('cH2', attachUnit(obj.cH2,'mol/m^3'));
    obj.m.param.set('cO2', attachUnit(obj.cO2,'mol/m^3'));
    obj.m.param.set('cFe', attachUnit(obj.cFe,'mol/m^3'));
    obj.m.param.set('cFe2p', attachUnit(obj.cFe2p,'mol/m^3'));

    % setup electrochemical surface reaction parameters
    for i=1:obj.nReactions 
        obj.m.param.set(obj.k0_id{i},attachUnit(obj.reactions{i}.k0,'m/s'),['standard rate constant of ',obj.reactions{i}.reaction]);
        obj.m.param.set(obj.E0_id{i},attachUnit(obj.reactions{i}.E0,'V'),['reduction potential of ',obj.reactions{i}.reaction]);
        obj.m.param.set(obj.n_id{i},obj.reactions{i}.n,['number of electrons participating in ',obj.reactions{i}.reaction]);
        obj.m.param.set(obj.beta_id{i},obj.reactions{i}.beta,['symmetry coefficient of ',obj.reactions{i}.reaction]);
    end

    for i=1:obj.numberOfSpecies
        for j = 1:obj.nReactions
            obj.m.param.set(obj.nu_id{i,j},num2str(obj.reactions{j}.flux(i)), sprintf('Flux of species %d due to reaction %s',i, obj.reactions{j}.reaction));
        end
    end
    %% create global functions
    fprintf('  Modifying global functions...\n');

    % analytical solution 
    obj.phi_pb.set('funcname', 'phi_pb');
    obj.phi_pb.set('args', {'x' 'phi0'});
    obj.phi_pb.set('expr', phiPBFunction);
    obj.phi_pb.set('fununit', '1');
    obj.phi_pb.set('argunit', '1');
    obj.phi_pb.set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
    obj.phi_pb.label('Potential distribution in symmetric electrolyte due to Poisson-Boltzmann model');

    obj.phi_pbx.set('funcname', 'phi_pbx');
    obj.phi_pbx.set('args', {'x' 'phi0'});
    %phi_pbx.set('expr', '- 4*RT/(abs(z1)*F_const)* tanh( abs(z1) * F_const* phi_bpe/(4*RT) )* exp(-x/lambdaD) / (lambdaD * (1- exp(-2*x/lambdaD) * tanh( abs(z1) * F_const * phi_bpe / (4*RT) )^2 ) )');
    obj.phi_pbx.set('expr', phiPBxFunction);
    obj.phi_pbx.set('fununit', '1');
    obj.phi_pbx.set('argunit', '1');
    obj.phi_pbx.set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
    obj.phi_pbx.label('Derivative of potential distribution');

    for i = 1:obj.numberOfSpecies   
        obj.c_pb{i}.set('funcname', obj.c_pb_id(i,:));
        obj.c_pb{i}.set('args', {'x' 'phi0'});
        obj.c_pb{i}.set('expr', cPBFunction{i});
        obj.c_pb{i}.set('fununit', '1');
        obj.c_pb{i}.set('argunit', '1');
        obj.c_pb{i}.set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
        %c_pb{i}.label('Concentration distribution in symmetric binary electrolyte due to Poisson-Boltzmann model');
    end

    %% modifying geometry
    fprintf('  Modifying geometry...\n');

    obj.geom.repairTol(1e-16); % should be smaller than smallest elements
    % obj.space.set('intervals','many');
    % space.set('p','0,lambdaD,L'); % an interval with three points
    % space.set('p','0,lambdaD/L,1'); % an interval with three points. dimensionless

    %electrode_surface = space.

    obj.space.set('size', {'W/L' 'epsilon'});
    obj.space.set('pos', {'-W/(2*L)' '0'});
%     obj.space.setIndex('layer', 'epsilon', 0);
%     obj.space.set('layername', {'dl'}); % diffuse layer
    obj.ddl.set('size', {'W/L' '1-epsilon'});
    obj.space.set('pos', {'-W/(2*L)' 'epsilon'});
    obj.leftEdgeOfBpe.setIndex('p', '-Wbpe/(2*L)', 0, 0);
    obj.leftEdgeOfBpe.setIndex('p', '0', 1, 0);
    obj.rightEdgeOfBpe.setIndex('p', 'Wbpe/(2*L)', 0, 0);
    obj.rightEdgeOfBpe.setIndex('p', '0', 1, 0);
    
    % obj.geom.feature('fin').set('repairtol', '1.0E-16');
    % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
    %obj.geom.feature('fin').set('repairtol', '1.0E-16');
	%obj.geom.run('fin');
    obj.geom.run;

    %% define selections
%     obj.surfaceNode.geom('geom',0); % point
%     obj.surfaceNode.set(1); % 1st created point as bpe surface
%     obj.surfaceNode.label('Surface node');
% 
%     obj.zetaNode.geom('geom',0); % point
%     obj.zetaNode.set(2); % 2nd created point as zeta plane
%     obj.zetaNode.label('Zeta plane node');
% 
%     obj.bulkBoundary.geom('geom',0); % point
%     obj.bulkBoundary.set(3); % 2rd created point far away at bulk solution
%     obj.bulkBoundary.label('Node at bulk');
% 
%     obj.regionOfFirstDebyeLength.geom('geom',1); % interval
%     obj.regionOfFirstDebyeLength.set(1); % ddl region of one debye length at surface
%     obj.regionOfFirstDebyeLength.label('DDL region');

    % point selections
    
    obj.leftBoundaryOfSurface.set('entitydim', '0');
    obj.leftBoundaryOfSurface.label('leftBoundaryOfSurface');
    obj.leftBoundaryOfSurface.set('xmin', '-W/(2*L)');
    obj.leftBoundaryOfSurface.set('xmax', '-W/(2*L)');
    obj.leftBoundaryOfSurface.set('ymin', '0');
    obj.leftBoundaryOfSurface.set('ymax', '0');
    
    obj.rightBoundaryOfSurface.set('entitydim', '0');
    obj.rightBoundaryOfSurface.label('rightBoundaryOfSurface');
    obj.rightBoundaryOfSurface.set('xmin', 'W/(2*L)');
    obj.rightBoundaryOfSurface.set('xmax', 'W/(2*L)');
    obj.rightBoundaryOfSurface.set('ymin', '0');
    obj.rightBoundaryOfSurface.set('ymax', '0');
    
    obj.leftBoundaryOfZetaPlane.set('entitydim', '0');
    obj.leftBoundaryOfZetaPlane.label('leftBoundaryOfZetaPlane');
    obj.leftBoundaryOfZetaPlane.set('xmin', '-W/(2*L)');
    obj.leftBoundaryOfZetaPlane.set('xmax', '-W/(2*L)');
    obj.leftBoundaryOfZetaPlane.set('ymin', 'epsilon');
    obj.leftBoundaryOfZetaPlane.set('ymax', 'epsilon');
    
    obj.rightBoundaryOfZetaPlane.set('entitydim', '0');
    obj.rightBoundaryOfZetaPlane.label('rightBoundaryOfZetaPlane');
    obj.rightBoundaryOfZetaPlane.set('xmin', 'W/(2*L)');
    obj.rightBoundaryOfZetaPlane.set('xmax', 'W/(2*L)');
    obj.rightBoundaryOfZetaPlane.set('ymin', 'epsilon');
    obj.rightBoundaryOfZetaPlane.set('ymax', 'epsilon');
    
    % edge selections

    obj.bpeSurface.set('entitydim', '1');
    obj.bpeSurface.label('bpeSurface');
    obj.bpeSurface.set('xmin', '0');
    obj.bpeSurface.set('xmax', '0');
    obj.bpeSurface.set('ymin', '0');
    obj.bpeSurface.set('ymax', '0');
    
    obj.bulkBoundary.set('entitydim', '1');
    obj.bulkBoundary.label('bulkBoundary');
    obj.bulkBoundary.set('xmin', '0');
    obj.bulkBoundary.set('xmax', '0');
    obj.bulkBoundary.set('ymin', '1');
    obj.bulkBoundary.set('ymax', '1');
    
    obj.zetaPlane.set('entitydim', '1');
    obj.zetaPlane.label('zetaPlane');
    obj.zetaPlane.set('xmin', '0');
    obj.zetaPlane.set('xmax', '0');
    obj.zetaPlane.set('ymin', '1/2*epsilon');
    obj.zetaPlane.set('ymax', 'epsilon');
   
    obj.electrodesAtFirstDebyeLength.set('entitydim', '1');
    obj.electrodesAtFirstDebyeLength.label('electrodes at first debye length');
    obj.electrodesAtFirstDebyeLength.set('ymin', 'epsilon/2');
    obj.electrodesAtFirstDebyeLength.set('ymax', 'epsilon/2');
    
    obj.electrodesAtRemainingDomain.set('entitydim', '1');
    obj.electrodesAtRemainingDomain.label('electrodes at remianing domain');
    obj.electrodesAtRemainingDomain.set('ymin', '1/2');
    obj.electrodesAtRemainingDomain.set('ymax', '1/2');
    
    obj.electrodes.set('entitydim', '1');
    obj.electrodes.set('input', {'electrodesAtFirstDebyeLength' 'electrodesAtRemainingDomain'});
    obj.electrodes.label('electrodes');
    
    obj.workingElectrode.set('entitydim', '1');
    obj.workingElectrode.label('workingElectrode');
    obj.workingElectrode.set('condition', 'inside');
    obj.workingElectrode.set('xmin', '-Inf');
    obj.workingElectrode.set('xmax', '-W/(2*L)');
    obj.workingElectrode.set('ymin', '-Inf');
    obj.workingElectrode.set('ymax', 'Inf');
    
    obj.counterElectrode.set('entitydim', '1');
    obj.counterElectrode.label('counterElectrode');
    obj.counterElectrode.set('condition', 'inside');
    obj.counterElectrode.set('xmin', 'W/(2*L)');
    obj.counterElectrode.set('xmax', 'Inf');
    obj.counterElectrode.set('ymin', '-Inf');
    obj.counterElectrode.set('ymax', 'Inf');
    
    obj.insulatorAdjacentToBpe.set('entitydim', '1');
    obj.insulatorAdjacentToBpe.set('input', {'bpeSurface'});
    obj.insulatorAdjacentToBpe.label('insulatorAdjacentToBpe');
        
    obj.entireSurface.set('entitydim', '1');
    obj.entireSurface.set('input', {'insulatorAdjacentToBpe' 'bpeSurface'});
    obj.entireSurface.label('entireSurface');
    
    % domain selections
    obj.regionOfFirstDebyeLength.set('entitydim', '1');
    obj.regionOfFirstDebyeLength.set('outputdim', '2');
    obj.regionOfFirstDebyeLength.set('input', {'bpeSurface'});
    obj.regionOfFirstDebyeLength.label('DDL region');

    %% local definitions
    fprintf('  Setting up local definitions...\n');
    % % operators
    % projectZetaPotential.selection.geom('geom', 0);
    % projectZetaPotential.selection.named('zetaNode');
    % projectZetaPotential.label('projectZetaPotential');
    % projectZetaPotential.set('opname', 'projectZetaPotential');
    % projectZetaPotential.selection('dstvertex1').geom('geom', 0);
    % projectZetaPotential.selection('dstvertex1').set(1);
    % %projectZetaPotential.selection('dstvertex1').named('surfaceNode');
    % %projectZetaPotential.selection('dstvertex2').geom('geom', 0);
    % %projectZetaPotential.selection('dstvertex3').geom('geom', 0);
    % projectZetaPotential.selection('srcvertex1').geom('geom', 0);
    % projectZetaPotential.selection('srcvertex1').set(2);
    % %projectZetaPotential.selection('srcvertex1').named('zetaNode');
    % %projectZetaPotential.selection('srcvertex2').geom('geom', 0);
    % %projectZetaPotential.selection('srcvertex3').geom('geom', 0);

    % operators
    
    leftBoundaryOfSurfaceVertexIndex = obj.leftBoundaryOfSurface.entities(0);
    rightBoundaryOfSurfaceVertexIndex = obj.rightBoundaryOfSurface.entities(0);
    leftBoundaryOfZetaPlaneVertexIndex = obj.leftBoundaryOfZetaPlane.entities(0);
    rightBoundaryOfZetaPlaneVertexIndex = obj.rightBoundaryOfZetaPlane.entities(0);

    obj.projectReactionPlaneToSurface.selection.geom('geom', 0);
    obj.projectReactionPlaneToSurface.selection.named('zetaPlane');
    obj.projectReactionPlaneToSurface.label('projectReactionPlaneToSurface');
    obj.projectReactionPlaneToSurface.set('opname', 'projectReactionPlaneToSurface');
%     obj.projectReactionPlaneToSurface.selection('dstvertex1').geom('geom', 0);
    obj.projectReactionPlaneToSurface.selection('dstvertex1').set(leftBoundaryOfSurfaceVertexIndex);
    obj.projectReactionPlaneToSurface.selection('dstvertex2').set(rightBoundaryOfSurfaceVertexIndex);
%     obj.projectReactionPlaneToSurface.selection('srcvertex1').geom('geom', 0);
    obj.projectReactionPlaneToSurface.selection('srcvertex1').set(leftBoundaryOfZetaPlaneVertexIndex);
    obj.projectReactionPlaneToSurface.selection('srcvertex2').set(rightBoundaryOfZetaPlaneVertexIndex);

    
    % project currents at surface to bulk exit
%     obj.projectSurfaceToBulkExit.selection.geom('geom', 0);
%     obj.projectSurfaceToBulkExit.selection.named('bulkBoundary');
%     obj.projectSurfaceToBulkExit.label('projectSurfaceToBulkExit');
%     obj.projectSurfaceToBulkExit.set('opname', 'projectSurfaceToBulkExit');
%     obj.projectSurfaceToBulkExit.selection('dstvertex1').geom('geom', 0);
%     obj.projectSurfaceToBulkExit.selection('dstvertex1').set(3);
%     obj.projectSurfaceToBulkExit.selection('srcvertex1').geom('geom', 0);
%     obj.projectSurfaceToBulkExit.selection('srcvertex1').set(1);
%     
    % integration operator on bpe surface
    obj.integrateSurface.selection.named('bpeSurface');
    obj.integrateSurface.label('integrateSurface');
    obj.integrateSurface.set('opname', 'integrateSurface');

    % integrate exit towards bulk
    obj.integrateBulkExit.selection.named('bulkBoundary');
    obj.integrateBulkExit.label('integrateBulkBoundary');
    obj.integrateBulkExit.set('opname', 'integrateBulkBoundary');

    % maximum on whole computational domain
    obj.maximumOnDomain.set('opname', 'maximumOnDomain');
    obj.maximumOnDomain.selection.all;

    %% define variables
    % concentrations and log concentrations
    obj.domainVariables.model('comp1');
    % for i = 1:numberOfSpecies   
    %     domainVariables.set(c_id(i,:), ...
    %          strcat('exp(',cellstr(C_id(i,:)),')'), ...
    %         'concentrations expressed in terms of log concentrations');
    % end
    obj.domainVariables.label('Variables valid over entire domain');


    obj.surfaceVariables.model('comp1');
    obj.surfaceVariables.selection.geom('geom', 1);
    obj.surfaceVariables.selection.named('bpeSurface');
    obj.surfaceVariables.label('Variables valid on surface');

    % for Dirichlet BC in potential:
    % surfaceVariables.set('V', 'phi_bpe-projectZetaPotential(phi)', 'Metal - reaction plane potential difference');
    % for Robin BC:
    % surfaceVariables.set('V', '(phi_bpe-phi)*UT', 'Metal - reaction plane potential difference, with dimension');
    obj.surfaceVariables.set('phi_s','phi_bpe');
    obj.surfaceVariables.set('V', '(phi_s-projectReactionPlaneToSurface(phi))*UT', 'Metal - reaction plane potential difference, with dimension');

    % % cathodic currents
    obj.surfaceVariables.set('i_cathodic', cathodicCurrent, 'Cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.');

    % % anodic currents    
    obj.surfaceVariables.set('i_anodic', anodicCurrent, 'Anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.');

    % totalCurrent = strjoin({cathodicCurrent,anodicCurrent},'');
    obj.surfaceVariables.set('i_total', totalCurrent, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.');

    % QUESTION: should n be includet in exponent???
    for i = 1:obj.nReactions
        obj.surfaceVariables.set(obj.i_id{i}, reactionCurrent{i}, ['current due to surface reaction ', obj.reactions{i}.reaction]);
    end

    % species Flux
    for i = 1:obj.numberOfSpecies
        obj.surfaceVariables.set(obj.N_id(i,:), fluxTerm(i));
    end
    
    obj.electrodeVariables.model('comp1');
    obj.electrodeVariables.selection.geom('geom', 1);
    obj.electrodeVariables.selection.named('electrodes');
    obj.electrodeVariables.label('Variables valid on both working and counter electrode sides');

    obj.weVariables.model('comp1');
    obj.weVariables.selection.geom('geom', 1);
    obj.weVariables.selection.named('workingElectrode');
    obj.weVariables.label('Variables valid  working electrode sides');

    obj.ceVariables.model('comp1');
    obj.ceVariables.selection.geom('geom', 1);
    obj.ceVariables.selection.named('counterElectrode');
    obj.ceVariables.label('Variables valid on counter electrode sides');
    
    obj.weVariables.set('phi_s', 'phi_we');
    obj.ceVariables.set('phi_s', 'phi_ce');

    
    %% Nernst Planck equations
    fprintf('  Setting up Nernst-Planck system...\n');

    % settings independent from transformation
    obj.NernstPlanckEquation.label('Nernst-Planck equation');
    obj.NernstPlanckEquation.prop('Units').set('DependentVariableQuantity', 'none');
    obj.NernstPlanckEquation.prop('Units').set('CustomSourceTermUnit', '1'); % really? check this

    obj.FluxAtSurface.selection.named('bpeSurface');
    obj.FluxAtSurface.label('Flux at surface of BPE');
    obj.FluxAtBulkNode.selection.named('bulkBoundary');
    obj.FluxAtBulkNode.label('Flux from bulk');

    obj.BulkConcentrationsBC.selection.named('bulkBoundary');
    obj.BulkConcentrationsBC.label('Fix bulk concentrations');

    switch obj.transformation
        case 'logConcentration'
            fprintf('    Nernst-Planck equation is transformed using loagrithmic concentrations.\n');
  
            % assume domainVariables already created
            for i = 1:numberOfSpecies   
                obj.domainVariables.set(obj.c_id{i}, ...
                     strcat('exp(',cellstr(obj.C_id{i}),')'), ...
                    'concentrations expressed in terms of log concentrations');
            end

            % add species concentration variables c
            obj.NernstPlanckEquation.field('dimensionless').component(obj.C_id);

            % set up units for variables (concentrations [c] = mol/m^3) and source term
            % (flux [R] = mol / (m^3 * s) )
            obj.NernstPlanckEquation.prop('Units').set('CustomDependentVariableUnit', '1');

            % define governing pde through diffusivities D and preset convection
            % coefficients determined by species charge and electrical moblity mu
            % NernstPlanckEquation.feature('gfeq1').set('Ga', prepArrayFromString( pdeGammaTerm,'%i',numberOfSpecies) );
            obj.NernstPlanckEquation.feature('gfeq1').set('Ga', pdeGammaTerm );

            %NernstPlanckEquation.feature('gfeq1').set('al', prepDiagFromString( pdeConvectionCoefficientAlpha, '%i', numberOfSpecies ));

            % sources in domain
            obj.NernstPlanckEquation.feature('gfeq1').set('f', pdeSourceTerm );

            % Nernst-Planck equation boundary conditions
            if( obj.bpeFluxOn )
                % set species flux at bpe surface due to surface reactions
                fprintf('    Surface species flux is activated.\n');

                obj.modifiedFluxTerm = strcat('surfaceFluxRampFactor*(',obj.N_id,')/exp(',obj.C_id,')');
                obj.FluxAtSurface.set('g', modifiedFluxTerm);
                obj.FluxAtBulkNode.set('g', cellstr(strcat('-integrateSurface(',modifiedFluxTerm,')')));
                % inflow at surface must correspond to outflow at bulk, and vice versa
            end

            % fix bulk concentration
            obj.BulkConcentrationsBC.set('r',strcat('log(',cellstr(obj.c_bulk_id),')') );

            obj.NernstPlanckEquation.feature('init1').set(obj.C_id(i,:), stract('log(',obj.c_pb_id{i},'(x))') );
        otherwise
            fprintf('    Nernst-Planck equation is used in its classical form.\n');

            % add species concentration variables c
            obj.NernstPlanckEquation.field('dimensionless').component(obj.c_id);

            % set up units for variables (concentrations [c] = mol/m^3) and source term
            % (flux [R] = mol / (m^3 * s) )
            obj.NernstPlanckEquation.prop('Units').set('CustomDependentVariableUnit', '1');

            % define governing pde through diffusivities D and preset convection
            % coefficients determined by species charge and electrical moblity mu
            obj.NernstPlanckEquation.feature('cfeq1').set('c',pdeNPDiffusionCoefficientC); % needs row array

            % OUTSOURCE THIS UGLY STUFF TO ANOTHER FILE
    %         pattern = diag(1:numberOfSpecies);
    %         pdeConvectionCoefficientAlphaMatrix = arrayfun( ... 
    %             @(p) iif(  p == 0, @() '0', true, @() pdeNPConvectionCoefficientAlpha{p} ), ...
    %             pattern, 'UniformOutput',false);
    %        temp = strdiag(pdeNPConvectionCoefficientAlpha);
    %         pdeConvectionCoefficientAlphaMatrix = reshape(pdeConvectionCoefficientAlphaMatrix,1,numberOfSpecies^2);
            pdeNPConvectionCoefficientAlphaXMatrix = flatten(strdiag(pdeNPConvectionCoefficientAlphaX));
            pdeNPConvectionCoefficientAlphaYMatrix = flatten(strdiag(pdeNPConvectionCoefficientAlphaY));

            pdeNPConvectionCoefficientAlphaMatrix = [pdeNPConvectionCoefficientAlphaXMatrix,pdeNPConvectionCoefficientAlphaYMatrix];
            %NernstPlanckEquation.feature('cfeq1').set('al', prepDiagFromString( pdeConvectionCoefficientAlpha, '%i', numberOfSpecies ));
            obj.NernstPlanckEquation.feature('cfeq1').set('al', pdeNPConvectionCoefficientAlphaMatrix);
            % 2d modification:
            % obj.NernstPlanckEquation.feature('cfeq1').set('al', {'z1*phix' 'z1*phiy'; '0' '0'; '0' '0'; 'z2*phix' 'z2*phiy'});


            % no sources in domain
            obj.NernstPlanckEquation.feature('cfeq1').set('f', pdeNPSourceTerm);

            % Nernst-Planck equation boundary conditions
            % set species flux at bulk to zero
            % FluxFromBulk = NernstPlanckEquation.create('FluxFromBulk', 'FluxBoundary', 0);
            % FluxFromBulk.selection.named('bulkBoundary');
            % FluxFromBulk.label('Flux from bulk');

            % set species flux at bpe surface due to surface reactions
            if( obj.bpeFluxOn )
                fprintf('    Surface species flux is activated.\n');
                % modifiedFluxTerm = strcat('surfaceFluxRampFactor*(',N_id,')');
                obj.FluxAtSurface.set('g', pdeNPbpeNeumannBC);
                % FluxAtBulkNode.set('g', cellstr(strcat('-integrateSurface(',pdeNPbpeNeumannBC,')'))); % necessary?

                % inflow at surface must correspond to outflow at bulk, and vice versa
                %  FluxAtSurface.set('g', cellstr(N_id));
                %  FluxAtBulkNode.set('g', cellstr(strcat('-integrateSurface(',N_id,')')));
                % inflow at surface must correspond to outflow at bulk, and vice versa

                %FluxAtSurface.set('g', zeros(numberOfSpecies,1)); % deactivate surface flux
            else
               obj.FluxAtSurface.set('g', zeros(obj.numberOfSpecies,1));
            end


            % fix bulk concentration
            obj.BulkConcentrationsBC.set('r', pdeNPbulkDirichletBC);

            % define initial values
            % modification depending on whether or not Stern layer thickness is part 
            % of computation domain
             for i = 1:obj.numberOfSpecies
                obj.NernstPlanckEquation.feature('init1').set(obj.c_id{i}, pdeNPinitialValues{i});
             end
    end

    %% Poisson equation
    fprintf('  Setting up Poisson equation...\n');

    obj.PoissonEquation.label('Poisson equation');

    % set single variable phi (electric potential)
    obj.PoissonEquation.field('dimensionless').component({'phi'});

    % set up units of variables ([phi] = V) and source term ([rho] = C/m^3)
    obj.PoissonEquation.prop('Units').set('DependentVariableQuantity', 'none');
    obj.PoissonEquation.prop('Units').set('CustomDependentVariableUnit', '1');
    obj.PoissonEquation.prop('Units').set('CustomSourceTermUnit', '1');

    % define governing pde through permeability
    obj.PoissonEquation.feature('cfeq1').set('c', pdePDiffusionCoefficientC);

     % prepare charge density term rho = sum rho_i = F sum Z_i c_i 
    % chargeDensityTerm = [ 'F_const*(', strjoin( cellstr( rcat( ...
    %     prepIndices('z',numberOfSpecies),...
    %     repmat('*',numberOfSpecies,1),...
    %     prepIndices('c',numberOfSpecies) ) )','+' ), ')'];
    obj.PoissonEquation.feature('cfeq1').set('f', pdePSourceTerm);

    % Poisson equation boundary conditions

    % set surface charge at boundaries, Robin BC
    % phi + lambdaS * phi' = phi_bpe, or
    % C_Stern*(phi_bpe - phi) = eps0*eps* phi' = - sigma
    % with lambdaS width of Stern layer, sigma surface charge density
    obj.bpeChargeDensityBC.selection.named('bpeSurface');
    obj.bpeChargeDensityBC.set('g', pdePbpeRobinBCg);
    obj.bpeChargeDensityBC.set('q', pdePbpeRobinBCq);
    obj.bpeChargeDensityBC.label('Charge density at surface');

    % potential at surface of bpe
    obj.bpePotentialBC.selection.named('bpeSurface');
    obj.bpePotentialBC.set('r', pdePbpeDirichletBC);
    obj.bpePotentialBC.label('potential at surface of BPE');

    %according initial values:
    obj.PoissonEquation.feature('init1').set('phi', pdePinitialValues);

    % chose, which bc to apply to Poisson equation at bpe surface
    switch obj.bpePoissonEquationBC 
        case 'Robin'
            fprintf('    Robin boundary conditions are used...\n');
            obj.bpeChargeDensityBC.active(true);
            obj.bpePotentialBC.active(false);

        case 'Dirichlet'
            fprintf('    Dirichlet boundary conditions are used...\n');
            obj.bpeChargeDensityBC.active(false);
            obj.bpePotentialBC.active(true);
    end

    % bulk net charge
    % bulkChargeDensityBC = PoissonEquation.create('bulkChargeDensityBC', 'FluxBoundary', 0);
    % bulkChargeDensityBC.selection.named('bulkBoundary');
    % bulkChargeDensityBC.set('g', '0');
    % bulkChargeDensityBC.label('Zero net charge density at bulk');

    % bulk potential
    obj.bulkPotentialBC.selection.named('bulkBoundary');
    obj.bulkPotentialBC.set('r', pdePbulkDirichletBC);
    obj.bulkPotentialBC.label('Bulk potential');
    
    % electrode potentials
    obj.electrodePotentialBC.selection.named('electrodes');
    obj.electrodePotentialBC.set('r', pdePelectrodeDirichletBC);
    obj.electrodePotentialBC.label('electrode potentials');
    % initial values, moved into switch fork above
    % PoissonEquation.feature('init1').set('phi', 'phi_bulk');

    %% Steady-state condition on BPE surface: zero net current
    obj.zeroSurfaceCurrent.label('zeroSurfaceCurrent');
    obj.zeroSurfaceCurrent.field('dimensionless').field('phi_bpe');
    obj.zeroSurfaceCurrent.field('dimensionless').component({'phi_bpe'});
    obj.zeroSurfaceCurrent.selection.named('bpeSurface');
    % zeroSurfaceCurrentEquation.selection.named('surfaceNode');

    % switch off standard equation:
    obj.zeroSurfaceCurrent.feature('dode1').set('f', '0'); % important, otherwise gives wrong results
    obj.zeroSurfaceCurrent.feature('dode1').set('da', '0');

    % not sure, which formulation is more exact
    %zeroSurfaceCurrentEquation.set('f', 'log(abs(i_anodic))-log(abs(i_cathodic))'); 
    obj.zeroSurfaceCurrentEquation.set('f', 'i_anodic+i_cathodic');
    obj.zeroSurfaceCurrentEquation.selection.named('bpeSurface');

    % initial value:
    obj.zeroSurfaceCurrent.feature('init1').set('phi_bpe', 'PHI_bpe/UT');

    %% meshing
    fprintf('  Setting mesh parameters...\n');

     % get mesh refinement from parameters
    obj.mesh1D();
    
    obj.standardMesh.label('standardMesh');
    %obj.standardMesh.feature('size').set('hauto', 2); % is this the automatic refinement?
   
    % mesh refinement at BPE surface, whole surface including insulated
    % pieces
    obj.bpeSurfaceMeshingEdge.selection.geom('geom',1);
    obj.bpeSurfaceMeshingEdge.selection.named('entireSurface');
    obj.bpeSurfaceMeshingEdge.label('Mesh refinement at BPE surface');

    obj.bpeSurfaceMeshingProperties.selection.geom('geom', 1); % what does this mean?
%   obj.bpeSurfaceMeshingProperties.selection.named('bpeSurface');
    obj.bpeSurfaceMeshingProperties.selection.named('entireSurface');
    obj.bpeSurfaceMeshingProperties.label('Mesh refinement at BPE surface');
    obj.bpeSurfaceMeshingProperties.set('custom', 'on');
    obj.bpeSurfaceMeshingProperties.set('hmaxactive', true);
    obj.bpeSurfaceMeshingProperties.set('hmax', min(obj.intFirstDebyeLength)); % maximum element width
    obj.bpeSurfaceMeshingProperties.set('hgradactive', true);
    obj.bpeSurfaceMeshingProperties.set('hgrad', min(obj.distributionFirstDebyeLength)); % maximum element growth coefficient

    % mesh size properties at zeta (or reaction) plane
    obj.zetaPlaneMeshingEdge.selection.geom('geom',1);
    obj.zetaPlaneMeshingEdge.selection.named('zetaPlane');
    obj.zetaPlaneMeshingEdge.label('Mesh refinement at zeta plane');
    
    obj.zetaPlaneMeshingProperties.selection.geom('geom', 1); % what does this mean?
%   obj.bpeSurfaceMeshingProperties.selection.named('bpeSurface');
    obj.zetaPlaneMeshingProperties.selection.named('zetaPlane');
    obj.zetaPlaneMeshingProperties.label('Mesh refinement at zeta plane');
    obj.zetaPlaneMeshingProperties.set('custom', 'on');
    obj.zetaPlaneMeshingProperties.set('hmaxactive', true);
    obj.zetaPlaneMeshingProperties.set('hmax', max(obj.intFirstDebyeLength)); % maximum element width
    obj.zetaPlaneMeshingProperties.set('hgradactive', true);
    obj.zetaPlaneMeshingProperties.set('hgrad', max(obj.distributionFirstDebyeLength)); % maximum element growth coefficient

    % mesh element distribution within diffuse layer
    obj.electrodeMeshingEdgeAtFirstDebyeLength.selection.geom('geom',1);
    obj.electrodeMeshingEdgeAtFirstDebyeLength.selection.named('electrodesAtFirstDebyeLength');
    obj.electrodeMeshingEdgeAtFirstDebyeLength.label('Mesh refinement at electrodes at first debye length above bpe');

    obj.electrodeMeshingPropertiesAtFirstDebyeLength.selection.geom('geom', 1); % what does this mean?
    obj.electrodeMeshingPropertiesAtFirstDebyeLength.selection.named('electrodesAtFirstDebyeLength');
    obj.electrodeMeshingPropertiesAtFirstDebyeLength.label('Mesh refinement at electrodes at first debye length above bpe');
    obj.electrodeMeshingPropertiesAtFirstDebyeLength.set('explicit', obj.distributionFirstDebyeLengthStr);
    obj.electrodeMeshingPropertiesAtFirstDebyeLength.set('type', 'explicit');

%     obj.electrodeMeshingPropertiesAtRemainingDomain.active(false);
%     obj.electrodeMeshingPropertiesAtRemainingDomain.selection.geom('geom', 1); % what does this mean?
%     obj.electrodeMeshingPropertiesAtRemainingDomain.selection.named('electrodesAtRemainingDomain');
%     obj.electrodeMeshingPropertiesAtRemainingDomain.label('Mesh refinement at electrodes at remaining geometry');
%     obj.electrodeMeshingPropertiesAtRemainingDomain.set('explicit', obj.distributionRemainingStr);
%     obj.electrodeMeshingPropertiesAtRemainingDomain.set('type', 'explicit');
%     
    % meshing properties for remaining geometry
%     obj.remainingGeometryMeshingProperties.label('Mesh proterties for remaining geometry');
%     obj.remainingGeometryMeshingProperties.set('custom', 'on');
%     obj.remainingGeometryMeshingProperties.set('hmaxactive', true);
%     obj.remainingGeometryMeshingProperties.set('hmax', '1/1000'); % 1000 seems to be a good value for water
%     obj.remainingGeometryMeshingProperties.set('hgradactive', true);

    % mesh size properties for all other boundaries, only set maximum size
    obj.remainingBoundariesMeshingEdge.selection.geom('geom',1);
    obj.remainingBoundariesMeshingEdge.selection.remaining();
    obj.remainingBoundariesMeshingEdge.label('Mesh refinement for remaining boundaries');
    
%     obj.remainingBoundariesMeshingProperties.selection.geom('geom', 1); % what does this mean?
%   obj.bpeSurfaceMeshingProperties.selection.named('bpeSurface');
%     obj.remainingBoundariesMeshingProperties.selection.remaining();
    obj.remainingBoundariesMeshingProperties.label('Mesh refinement for remaining boundaries');
    obj.remainingBoundariesMeshingProperties.set('custom', 'on');
    obj.remainingBoundariesMeshingProperties.set('hmaxactive', true);
    obj.remainingBoundariesMeshingProperties.set('hmax', max(obj.intRemaining)); % maximum element width
%     obj.zetaPlaneMeshingProperties.set('hgradactive', true);
%     obj.zetaPlaneMeshingProperties.set('hgrad', max(obj.distributionFirstDebyeLength)); % maximum element growth coefficient

    % mesh in remaining 2d domain
    obj.meshingDomain.selection.geom('geom', 2);
    obj.meshingDomain.selection.all;
    obj.meshingDomain.label('Free triangualar mesh');

    
    obj.domainMeshingProperties.selection.geom('geom', 2);
    obj.domainMeshingProperties.selection.all;
    obj.domainMeshingProperties.label('Free triangualar mesh maximum element size');
    obj.domainMeshingProperties.set('custom', 'on');
    obj.domainMeshingProperties.set('hmaxactive', true);
    obj.domainMeshingProperties.set('hmax', max(obj.intRemaining));
    
    % testing mesh
    obj.testingMesh.label('testingMesh');
    %obj.standardMesh.feature('size').set('hauto', 2); % is this the automatic refinement?
   
    % mesh refinement at BPE surface, whole surface including insulated
    % pieces
    obj.testingMeshBpeSurfaceMeshingEdge.selection.geom('geom',1);
    obj.testingMeshBpeSurfaceMeshingEdge.selection.named('entireSurface');
    obj.testingMeshBpeSurfaceMeshingEdge.label('Mesh refinement at BPE surface');

    obj.testingMeshBpeSurfaceMeshingProperties.selection.geom('geom', 1); % what does this mean?
    obj.testingMeshBpeSurfaceMeshingProperties.selection.named('entireSurface');
    obj.testingMeshBpeSurfaceMeshingProperties.label('Mesh refinement at BPE surface');
    obj.testingMeshBpeSurfaceMeshingProperties.set('custom', 'on');
    obj.testingMeshBpeSurfaceMeshingProperties.set('hmaxactive', true);
    obj.testingMeshBpeSurfaceMeshingProperties.set('hmax', 'epsilon/10'); % maximum element width
    obj.testingMeshBpeSurfaceMeshingProperties.set('hgradactive', true);
    obj.testingMeshBpeSurfaceMeshingProperties.set('hgrad', '1.01'); % maximum element growth coefficient

    % mesh in remaining 2d domain
    obj.testingMeshTriangular.selection.geom('geom', 2);
    obj.testingMeshTriangular.selection.all;
    obj.testingMeshTriangular.label('Free triangualar mesh for testing purposes');
    
    obj.testingMeshTriangularProperties.selection.geom('geom', 2);
    obj.testingMeshTriangularProperties.selection.all;
    obj.testingMeshTriangularProperties.label('Free triangualar mesh for testing purposes maximum element size');
    obj.testingMeshTriangularProperties.set('custom', 'on');
    obj.testingMeshTriangularProperties.set('hmaxactive', true);
    obj.testingMeshTriangularProperties.set('hmax', '1/4');
    % obj.standardMesh.run;
    
    %% update views
%     obj.updateViews();
    %% update datasets
%     obj.bpeSurfaceResults.label('bpeSurfaceResults');
%     obj.bpeSurfaceResults.selection.named('bpeSurface');
% 
%     obj.entireSurfaceResults.label('entireSurfaceResults');
%     obj.entireSurfaceResults.selection.named('entireSurface');
% 
%     obj.zetaPlaneResults.label('zetaPlaneResults');
%     obj.zetaPlaneResults.selection.named('zetaPlane');
% 
%     obj.bulkBoundaryResults.label('bulkBoundaryResults');
%     obj.bulkBoundaryResults.selection.named('bulkBoundary');
% 
%     obj.weResults.label('weResults');
%     obj.weResults.selection.named('workingElectrode');
% 
%     obj.ceResults.label('ceResults');
%     obj.ceResults.selection.named('counterElectrode');
% 
%     obj.centralCrossectionResults.label('centralCrossectionRseults');
%     obj.centralCrossectionResults .set('genpoints', {'0' '0'; '0' '1'});
% 
%     obj.centralDDLCrossectionResults.label('centralDDLCrossectionResults');
%     obj.centralDDLCrossectionResults.set('genpoints', {'0' '0'; '0' num2str(obj.epsilon)});
end