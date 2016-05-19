function obj = updatePhysicsForWeakForm(obj)
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
        if ~isempty(fluxDirections)
            nonZeroFluxSummand = specificFluxSummand(i,fluxDirections(i,:)~=0);
            fluxTerm{i} = strjoin( nonZeroFluxSummand, '' );
            if strcmp(fluxTerm{i},'')
                fluxTerm{i} = '0';
            end
        else
            fluxTerm{i} = '0';
        end
    end
    
     % dimensionless fluxes
    fluxTermDimless = prepTerm('L/(D_id*c_ref)*N_id','D_id','N_id',obj.D_id,obj.N_id);

    % fluxTerm{ strcmp(fluxTerm,'') } = '0';

    % concentrations of species participating in surface reactions
%     oxidantConcentrations = cellfun( @(r) r.oxidantConcentrations{1}, obj.reactions','UniformOutput',false);
%     reductantConcentrations = cellfun( @(r) r.reductantConcentrations{1}, obj.reactions','UniformOutput',false);

    % IMPORTANT 
    % opposite of convention in comsol
    %   anodic currents NEGATIVE, cathodic currents POSITIVE

    % comsol convention:
    %   anodic currents POSITIVE, cathodic currents NEGATIVE

    % anodic currents 
    
%     anodicTerms = prepTerm(...
%         'n_id*F_const*k0_id*(c_R*exp(n_id*(1-beta_id)*F_const/RT*(V-E0_id)))',...
%         'n_id','k0_id','c_R','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,reductantConcentrations,obj.beta_id,obj.E0_id);
    anodicTerms = cellfun( @(r) r.anodicCurrentTerm, obj.reactions, 'UniformOutput', false);
    emptyAnodicTerms = cellfun( @(c) isempty(c), anodicTerms );
%     anodicTerms{emptyTerms} = {};
    anodicCurrent = strjoin(anodicTerms(~emptyAnodicTerms),'+');

    % cathodic currents
%     cathodicTerms = prepTerm(...
%         '-n_id*F_const*k0_id*(c_Ox*exp(-n_id*beta_id*F_const/RT*(V-E0_id)))',...
%         'n_id','k0_id','c_Ox','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,oxidantConcentrations,obj.beta_id,obj.E0_id);    
    cathodicTerms = cellfun( @(r) r.cathodicCurrentTerm, obj.reactions, 'UniformOutput', false );
    emptyCathodicTerms = cellfun( @(c) isempty(c), cathodicTerms );
%     cathodicTerms(emptyTerms) = '';
    cathodicCurrent = strjoin(cathodicTerms(~emptyCathodicTerms),'');

    totalCurrent = strjoin({anodicCurrent,cathodicCurrent},'');

    % current per reaction
    % Butler-Volmer, cathodic direction > 0 in form of, as in Bard p. 96
    % i = n*F/RT*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
    % in anodic direction > 0:
%     reactionCurrent = prepTerm(...
%         'n_id*F_const*k0_id*( c_R*exp(n_id*(1-beta_id)*F_const/RT*(V-E0_id)) - c_Ox*exp(-n_id*beta_id*F_const/RT*(V-E0_id)) ) ',...
%         'n_id','k0_id','c_Ox','c_R','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,oxidantConcentrations,reductantConcentrations,obj.beta_id,obj.E0_id);
    anodicTerms{emptyAnodicTerms} = '';
    cathodicTerms{emptyCathodicTerms} = '';
    reactionCurrent = prepTerm('AnodicCathodic','Anodic','Cathodic',anodicTerms',cathodicTerms');
    %         n_id{i},'*F_const*',k0_id{i},'*( ',...
    %         reactions{i}.oxidantConcentrations{1},...
    %         '*exp(-',beta_id{i},'*F_const/RT * (V-',E0_id{i},'))',...
    %         '-',...
    %         reactions{i}.reductantConcentrations{1},...
    %         '*exp((1-',beta_id{i},')*F_const/RT * (V-',E0_id{i},')) )'...
    %         ), ['current due to surface reaction ', reactions{i}.reaction]);
 
    % dimensionless currents
    reactionCurrentDimless  = prepTerm('L*i_id/(D_ref*c_ref*F_const)','i_id',obj.i_id);
    anodicCurrentDimless    = 'L*i_anodic/(D_ref*c_ref*F_const)';
    cathodicCurrentDimless  = 'L*i_cathodic/(D_ref*c_ref*F_const)';
    totalCurrentDimless     = 'L*i_total/(D_ref*c_ref*F_const)';


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
%     pdeNPbulkDirichletBC    = prepTerm('c_bulk_id/c_ref','c_bulk_id',obj.c_bulk_id); % fixed concentration at bulk
    pdeNPbulkDirichletBC    = prepTerm('c_0_id','c_0_id',obj.c_0_id); % fixed concentration at bulk
%     pdeNPbpeNeumannBC       = prepTerm('L*N_id/(D_id*c_ref)','N_id','D_id',obj.N_id,obj.D_id); % species flux at surface
    pdeNPbpeNeumannBC       = prepTerm('smoothenBpeBC(x)*N_dimless_id','N_dimless_id',obj.N_dimless_id); % species flux at surface
    pdeNPbpeNeumannBC1d     = prepTerm('N_dimless_id','N_dimless_id',obj.N_dimless_id); % species flux at surface


    % initial value terms
%     pdeNPinitialValues = prepTerm('c_pb_id(y,phi_bpe)','c_pb_id',obj.c_pb_id);
    pdeNPinitialValues = prepTerm('c_0_id','c_0_id',obj.c_0_id);

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

   
    %phi0_term = 'PHI_bpe/UT'; % dimensionless surface potential

    % pdePbulkDirichletBC   = 'phi_bulk/UT'; % voltage at bulk, usually zero
%     pdePbulkDirichletBC   = 'phi_we-x*L/W*deltaPhi'; % voltage at bulk, usually zero
%     pdePbulkDirichletBC   = 'phi_pb(1,phi_we)-x*L/W*(phi_pb(1,phi_we)-phi_pb(1,phi_ce))'; % voltage at bulk, usually zero
%     pdePbulkDirichletBC = '0';
    pdePbulkDirichletBC = 'phi_s';
    
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
   
    pdePbpeRobinBCg1d     = 'C_Stern*phi_s';
    pdePbpeRobinBCq1d     = 'C_Stern';
    pdePbpeRobinBCg     = 'smoothenBpeBC(x)*C_Stern*phi_s';
    pdePbpeRobinBCq     = 'smoothenBpeBC(x)*C_Stern';

    % for delta = 0 => lambdaS = 0, BC reduces to
    %   phi = phi0

    pdePbpeDirichletBC1d = '*phi_s';
    pdePbpeDirichletBC = 'smoothenBpeBC(x)*phi_s';

    % for working and counter electrodes
%     pdePelectrodeDirichletBC = 'phi_pb(y,phi_s)';
    pdePelectrodeDirichletBC = 'phi_s';

    % initial values
    % pdePinitialValues = 'phi_pb(y,phi_bpe)';
    pdePinitialValues = 'phi0';



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
  
    % dimensionless flux N*:
    %   N   = - (D grad c + z F u c grad phi)
    %       = - c_ref D/L (grad* c* + z c* grad* phi*)
    %   N L / (c_ref D_ref) = - D / D_ref (grad* c* + z c* grad* phi*) = N*
    
    nxTerm = prepTerm('(-cx_id-z_id*c_id*phix)','D_id','cx_id','z_id','c_id',obj.D_id, obj.cx_id,obj.z_id,obj.c_id);
    nyTerm = prepTerm('(-cy_id-z_id*c_id*phiy)','D_id','cy_id','z_id','c_id',obj.D_id, obj.cy_id,obj.z_id,obj.c_id);
    
    NxTerm = prepTerm('(D_id*c_ref/L*nx_id)','D_id','nx_id',obj.D_id,obj.nx_id); % dimensional
    NyTerm = prepTerm('(D_id*c_ref/L*ny_id)','D_id','ny_id',obj.D_id,obj.ny_id);

    ixSummand = prepTerm('z_id*nx_id','z_id','nx_id',obj.z_id,obj.nx_id);
    iySummand = prepTerm('z_id*ny_id','z_id','ny_id',obj.z_id,obj.ny_id);
    
    IxTerm = 'D_ref*c_ref*F_const/L*ix';
    IyTerm = 'D_ref*c_ref*F_const/L*iy';
    
    ixTerm = strcat(strjoin(ixSummand','+' ));
    iyTerm = strcat(strjoin(iySummand','+' ));
    
    % local ionic conductivity
    kappaSummand = prepTerm('z_id^2*c_id*D_id','z_id','c_id','D_id',obj.z_id,obj.c_id,obj.D_id);
    kappaTerm    = ['F_const/RT*(',strcat(strjoin(kappaSummand','+' )),')'];

    c0Term = prepTerm('c_bulk_id/c_ref','c_bulk_id',obj.c_bulk_id);
%     phi0Term = '0';
    phi0Term = '-deltaPhi*x/w';
    
    %% 0d component
    obj.m.variable('var0d').set('i_total', totalCurrent);
%     obj.m.variable('var0d').set('phi_bpe', 'V');
    obj.m.physics('zeroNetCurrent0d').feature('ge1').setIndex('name', 'V', 0, 0);
	obj.m.physics('zeroNetCurrent0d').feature('ge1').setIndex('equation', 'i_total', 0, 0);
    
    %% 1d component
    c0Term1d = prepTerm('c_pb_id(x,phi_bpe_init)','c_pb_id',obj.c_pb_id);
   
% weak formulation
    wNernstPlanck1d               = prepTerm('(-c_idx-z_id*c_id*phix)*test(c_idx)','c_id','z_id',obj.c_id,obj.z_id);
    wPoisson1d                    = prepTerm('-epsilon^2*phix*test(phix)+chargeDensityRampFactor*chargeDensityTerm*test(phi)','chargeDensityTerm',chargeDensityTerm);

    wSurfaceFluxBC1d              = prepTerm('surfaceFluxRampFactor*NN_id*test(c_id)','NN_id','c_id',obj.N_dimless_id,obj.c_id);
    wBulkFluxBC1d                 = prepTerm('surfaceFluxRampFactor*integrateSurface1d(NN_id)*test(c_id)','NN_id','c_id',obj.N_dimless_id,obj.c_id);
    wSurfaceChargeDensityBC1d     = '-chargeDensityRampFactor*epsilon/delta*(phi-phi_s)*test(phi)';


    wBulkConcentrationBC1d        = prepTerm('-lambdaCBulk_id*test(c_id)-test(lambdaCBulk_id)*(c_id-c_id_0)','lambdaCBulk_id','c_id',obj.lambdaCBulk_id,obj.c_id);
    wBulkPotentialBC1d            = '-lambdaPhiBulk*test(phi)-test(lambdaPhiBulk)*phi';
    %wElectrodePotentialBC1d       = '-lambdaPhiBulk*test(phi)-test(lambdaPhiBulk)*(phi-phi_s)';

    wEquations1d = [wNernstPlanck1d; wPoisson1d];

% variables
    
    obj.m.variable('DomainVar1d').label('DomainVar1d');
    obj.m.variable('DomainVar1d').model('comp1d');
    obj.m.variable('DomainVar1d').selection.geom('geom1d', 1);
    obj.m.variable('DomainVar1d').selection.all;
      
    obj.m.variable('DomainVar1d').set('phi0', 'phi_pb(x,phi_bpe_init)');
    
    for i = 1:obj.numberOfSpecies
        obj.m.variable('DomainVar1d').set(obj.c_0_id{i}, c0Term1d{i});
%         obj.m.variable('BulkVar1d').set(obj.c_0_id{i}, c0Term1d{i});
    end
    
    for i = 1:obj.numberOfSpecies
       obj.m.variable('DomainVar1d').set(obj.nx_id{i},nxTerm{i});
       obj.m.variable('DomainVar1d').set(obj.Nx_id{i},NxTerm{i});
    end
    
    obj.m.variable('DomainVar1d').set('ix',ixTerm);
    obj.m.variable('DomainVar1d').set('Ix',IxTerm);
    obj.m.variable('DomainVar1d').set('kappa_loc',kappaTerm);

    % surface variable
    obj.m.variable('SurfaceVar1d').model('comp1d');
    obj.m.variable('SurfaceVar1d').selection.geom('geom1d', 0);
    obj.m.variable('SurfaceVar1d').selection.named('surfaceVertex1d');
    obj.m.variable('SurfaceVar1d').label('SurfaceVar1d');
    
    obj.m.variable('SurfaceVar1d').set('phi_s','phi_bpe');
    obj.m.variable('SurfaceVar1d').set('V', '(phi_s-phi)*UT', 'Metal - reaction plane potential difference, with dimension');
    
    obj.m.variable('SurfaceVar1d').set('i_cathodic', cathodicCurrent, 'Cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.');
    obj.m.variable('SurfaceVar1d').set('ii_cathodic', cathodicCurrentDimless, 'Dimensionless cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.');

    % % anodic currents    
    obj.m.variable('SurfaceVar1d').set('i_anodic', anodicCurrent, 'Anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.');
    obj.m.variable('SurfaceVar1d').set('ii_anodic', anodicCurrentDimless, 'Dimensionless anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.');

    % totalCurrent = strjoin({cathodicCurrent,anodicCurrent},'');
    obj.m.variable('SurfaceVar1d').set('i_total', totalCurrent, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.');
    obj.m.variable('SurfaceVar1d').set('ii_total', totalCurrentDimless, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.');

    for i = 1:obj.nReactions
        obj.m.variable('SurfaceVar1d').set(obj.i_id{i}, reactionCurrent{i}, ['current due to surface reaction ', obj.reactions{i}.reaction]);
        obj.m.variable('SurfaceVar1d').set(obj.i_dimless_id{i}, reactionCurrentDimless{i}, ['dimensionless current due to surface reaction ', obj.reactions{i}.reaction]);
    end

    % species Flux
    for i = 1:obj.numberOfSpecies
        obj.m.variable('SurfaceVar1d').set(obj.N_id{i}, fluxTerm{i});
        obj.m.variable('SurfaceVar1d').set(obj.N_dimless_id{i}, fluxTermDimless{i});
    end

    obj.m.variable('BulkVar1d').model('comp1d');
    obj.m.variable('BulkVar1d').selection.geom('geom1d', 0);
    obj.m.variable('BulkVar1d').selection.named('bulkExitVertex1d');
    obj.m.variable('BulkVar1d').label('BulkVar1d');
    obj.m.variable('BulkVar1d').set('phi_s','phi_bulk');

    
    % Nernst Planck equations
    
% weak formulation, 1d
    % set discretization
    obj.m.physics('WeakFormulation1d').prop('ShapeProperty').set('order', '4');
    obj.m.physics('WeakFormulation1d').prop('ShapeProperty').set('valueType', 'real');
    obj.m.physics('WeakFormulation1d').prop('ShapeProperty').set('shapeFunctionType', 'shherm');
    
    obj.m.physics('WeakFormulation1d').field('dimensionless').component([obj.c_id {'phi'}]);
    obj.m.physics('WeakFormulation1d').feature('wfeq1').set('weak', wEquations1d);
  
    % Nernst Planck
    for i=1:obj.numberOfSpecies
        % initial values
        obj.m.physics('WeakFormulation1d').feature('init1').set(obj.c_id{i},obj.c_0_id{i});
        
        % surface flux
        obj.m.physics('WeakFormulation1d').feature(obj.surfaceFluxBC_id{i}).selection.named('surfaceVertex1d');
        obj.m.physics('WeakFormulation1d').feature(obj.surfaceFluxBC_id{i}).set('weakExpression',wSurfaceFluxBC1d{i});
        
        % bulk flux (stabilization)
        obj.m.physics('WeakFormulation1d').feature(obj.bulkFluxBC_id{i}).selection.named('bulkExitVertex1d');
        obj.m.physics('WeakFormulation1d').feature(obj.bulkFluxBC_id{i}).set('weakExpression',wBulkFluxBC1d{i});
        
        % bulk concentration
        obj.m.physics('WeakFormulation1d').feature(obj.bulkConcentrationBC_id{i}).selection.named('bulkExitVertex1d');
        obj.m.physics('WeakFormulation1d').feature(obj.bulkConcentrationBC_id{i}).set('weakExpression',wBulkConcentrationBC1d{i});
        obj.m.physics('WeakFormulation1d').feature(obj.bulkConcentrationBC_id{i}).feature([obj.bulkConcentrationBC_id{i},'_aux']).set('fieldVariableName',obj.lambdaCBulk_id{i});
    end
       
    % Poisson
    
    % initial value
    obj.m.physics('WeakFormulation1d').feature('init1').set('phi', 'phi0');
    
    % bulk potential
    obj.m.physics('WeakFormulation1d').feature('BulkPotentialBC').selection.named('bulkExitVertex1d');
    obj.m.physics('WeakFormulation1d').feature('BulkPotentialBC').set('weakExpression',wBulkPotentialBC1d);
    obj.m.physics('WeakFormulation1d').feature('BulkPotentialBC').feature('BulkPotentialBC_aux').set('fieldVariableName','lambdaPhiBulk');
    
    % surface charge density
    obj.m.physics('WeakFormulation1d').feature('SurfaceChargeDensityBC').selection.named('surfaceVertex1d');
    obj.m.physics('WeakFormulation1d').feature('SurfaceChargeDensityBC').set('weakExpression',wSurfaceChargeDensityBC1d);


% classical PDE formulation

%     % settings independent from transformation
%     obj.m.physics('NernstPlanckEquation1d').label('NernstPlanckEquation1d');
%     obj.m.physics('NernstPlanckEquation1d').prop('Units').set('DependentVariableQuantity', 'none');
%     obj.m.physics('NernstPlanckEquation1d').prop('Units').set('CustomSourceTermUnit', '1'); % really? check this
% 
%     obj.m.physics('NernstPlanckEquation1d').feature('FluxAtSurface1d').selection.named('surfaceVertex1d');
%     obj.m.physics('NernstPlanckEquation1d').feature('FluxAtSurface1d').label('FluxAtSurface1d');
% 
%     obj.m.physics('NernstPlanckEquation1d').feature('BulkConcentrationsBC1d').selection.named('bulkExitVertex1d');
%     obj.m.physics('NernstPlanckEquation1d').feature('BulkConcentrationsBC1d').label('BulkConcentrationsBC1d');
% 
%     % add species concentration variables c
%     obj.m.physics('NernstPlanckEquation1d').field('dimensionless').component(obj.c_id);
% 
%     % set up units for variables (concentrations [c] = mol/m^3) and source term
%     % (flux [R] = mol / (m^3 * s) )
%     obj.m.physics('NernstPlanckEquation1d').prop('Units').set('CustomDependentVariableUnit', '1');
% 
%     % define governing pde through diffusivities D and preset convection
%     % coefficients determined by species charge and electrical moblity mu
%     obj.m.physics('NernstPlanckEquation1d').feature('cfeq1').set('c',pdeNPDiffusionCoefficientC); % needs row array
% 
%     pdeNPConvectionCoefficientAlphaXMatrix = flatten(strdiag(pdeNPConvectionCoefficientAlphaX));
%     obj.m.physics('NernstPlanckEquation1d').feature('cfeq1').set('al', pdeNPConvectionCoefficientAlphaXMatrix);
%  
%     % no sources in domain
%     obj.m.physics('NernstPlanckEquation1d').feature('cfeq1').set('f', pdeNPSourceTerm);
% 
%     % species flux at surface
%     obj.m.physics('NernstPlanckEquation1d').feature('FluxAtSurface1d').set('g', pdeNPbpeNeumannBC1d);
% 
%     % fix bulk concentration
%     obj.m.physics('NernstPlanckEquation1d').feature('BulkConcentrationsBC1d').set('r', pdeNPbulkDirichletBC);
% 
%     % define initial values
%     % modification depending on whether or not Stern layer thickness is part 
%     % of computation domain
%      for i = 1:obj.numberOfSpecies
%         obj.m.physics('NernstPlanckEquation1d').feature('init1').set(obj.c_id{i}, pdeNPinitialValues{i});
%      end        
% 
%     % Poisson equation
%     obj.m.physics('PoissonEquation1d').label('PoissonEquation1d');
% 
%     % set single variable phi (electric potential)
%     obj.m.physics('PoissonEquation1d').field('dimensionless').component({'phi'});
% 
%     % set up units of variables ([phi] = V) and source term ([rho] = C/m^3)
%     obj.m.physics('PoissonEquation1d').prop('Units').set('DependentVariableQuantity', 'none');
%     obj.m.physics('PoissonEquation1d').prop('Units').set('CustomDependentVariableUnit', '1');
%     obj.m.physics('PoissonEquation1d').prop('Units').set('CustomSourceTermUnit', '1');
% 
%     % define governing pde through permeability
%     obj.m.physics('PoissonEquation1d').feature('cfeq1').set('c', pdePDiffusionCoefficientC);
%     obj.m.physics('PoissonEquation1d').feature('cfeq1').set('f', pdePSourceTerm);
% 
%     % Poisson equation boundary conditions
%     
%     % set surface charge at boundaries, Robin BC
%     % phi + lambdaS * phi' = phi_bpe, or
%     % C_Stern*(phi_bpe - phi) = eps0*eps* phi' = - sigma
%     % with lambdaS width of Stern layer, sigma surface charge density
%     obj.m.physics('PoissonEquation1d').feature('bpeChargeDensityBC1d').selection.named('surfaceVertex1d');
%     obj.m.physics('PoissonEquation1d').feature('bpeChargeDensityBC1d').set('g', pdePbpeRobinBCg1d);
%     obj.m.physics('PoissonEquation1d').feature('bpeChargeDensityBC1d').set('q', pdePbpeRobinBCq1d);
%     obj.m.physics('PoissonEquation1d').feature('bpeChargeDensityBC1d').label('bpeChargeDensityBC1d');
% 
%     % potential at surface of bpe
%     obj.m.physics('PoissonEquation1d').feature('bpePotentialBC1d').selection.named('surfaceVertex1d');
%     obj.m.physics('PoissonEquation1d').feature('bpePotentialBC1d').set('r', pdePbpeDirichletBC1d);
%     obj.m.physics('PoissonEquation1d').feature('bpePotentialBC1d').label('bpePotentialBC1d');
%     
% 
%     %according initial values:
%     obj.m.physics('PoissonEquation1d').feature('init1').set('phi', pdePinitialValues);
% 
%     % bulk potential
%     obj.m.physics('PoissonEquation1d').feature('bulkPotentialBC1d').selection.named('bulkExitVertex1d');
%     obj.m.physics('PoissonEquation1d').feature('bulkPotentialBC1d').set('r', pdePbulkDirichletBC);
%     obj.m.physics('PoissonEquation1d').feature('bulkPotentialBC1d').label('bulkPotentialBC1d');
%  

    % Steady-state condition on BPE surface: zero net current
    obj.m.physics('zeroNetCurrent1d').label('zeroNetCurrent1d');
    obj.m.physics('zeroNetCurrent1d').field('dimensionless').field('phi_bpe');
    obj.m.physics('zeroNetCurrent1d').field('dimensionless').component({'phi_bpe'});
    obj.m.physics('zeroNetCurrent1d').selection.named('surfaceVertex1d');
    % zeroSurfaceCurrentEquation.selection.named('surfaceNode');

    % switch off standard equation:
    obj.m.physics('zeroNetCurrent1d').feature('dode1').set('f', '0'); % important, otherwise gives wrong results
    obj.m.physics('zeroNetCurrent1d').feature('dode1').set('da', '0');

    % not sure, which formulation is more exact
    %zeroSurfaceCurrentEquation.set('f', 'log(abs(i_anodic))-log(abs(i_cathodic))'); 
    obj.m.physics('zeroNetCurrent1d').feature('zeroSurfaceCurrentEquation1d').set('f', 'i_anodic+i_cathodic');
    obj.m.physics('zeroNetCurrent1d').feature('zeroSurfaceCurrentEquation1d').selection.named('surfaceVertex1d');

    % initial value:
    obj.m.physics('zeroNetCurrent1d').feature('init1').set('phi_bpe', 'phi_bpe_init');
    
    %% 2d component, weak formulation
    wNernstPlanck               = prepTerm('(-c_idx-z_id*c_id*phix)*test(c_idx)+(-c_idy-z_id*c_id*phiy)*test(c_idy)','c_id','z_id',obj.c_id,obj.z_id);
    wPoisson                    = prepTerm('-epsilon^2*(phix*test(phix)+phiy*test(phiy))+chargeDensityRampFactor*chargeDensityTerm*test(phi)','chargeDensityTerm',chargeDensityTerm);
    
    wSurfaceFluxBC              = prepTerm('surfaceFluxRampFactor*smoothenBpeBC(x)*NN_id*test(c_id)','NN_id','c_id',obj.N_dimless_id,obj.c_id);
    %wNPBuFluxBC     = prepTerm('NN_id*test(c_id)',{'NN_id','c_id'},{obj.NN_id,obj.c_id});
    wSurfaceChargeDensityBC     = '-chargeDensityRampFactor*epsilon/delta*(phi-phi_s)*test(phi)';
    
    
    wBulkConcentrationBC        = prepTerm('-lambdaCBulk_id*test(c_id)-test(lambdaCBulk_id)*(c_id-c_id_0)','lambdaCBulk_id','c_id',obj.lambdaCBulk_id,obj.c_id);
%     wBulkPotentialBC            = '-lambdaPhiBulkExit*test(phi)-test(lambdaPhiBulkExit)*phi';
    wElectrodePotentialBC       = '-lambdaPhiElectrode*test(phi)-test(lambdaPhiElectrode)*(phi-phi_s)';

    wEquations = [wNernstPlanck; wPoisson];
    
% define variables 2d
    
    % domain variables
    
    % concentrations and log concentrations
    obj.domainVariables.model('comp1');
    obj.domainVariables.selection.geom('geom', 2);
    obj.domainVariables.selection.all;
    % for i = 1:numberOfSpecies   
    %     domainVariables.set(c_id(i,:), ...
    %          strcat('exp(',cellstr(C_id(i,:)),')'), ...
    %         'concentrations expressed in terms of log concentrations');
    % end
    
    obj.domainVariables.set('phi0', phi0Term);
    for i = 1:obj.numberOfSpecies
        obj.domainVariables.set(obj.c_0_id{i}, c0Term{i});
    end
    
    obj.domainVariables.label('Variables valid over entire domain');
    
    for i = 1:obj.numberOfSpecies
            obj.domainVariables.set(obj.nx_id{i},nxTerm{i});
            obj.domainVariables.set(obj.ny_id{i},nyTerm{i});
            obj.domainVariables.set(obj.Nx_id{i},NxTerm{i});
            obj.domainVariables.set(obj.Ny_id{i},NyTerm{i});
    end
    obj.domainVariables.set('ix',ixTerm);
    obj.domainVariables.set('iy',iyTerm);
    obj.domainVariables.set('Ix',IxTerm);
    obj.domainVariables.set('Iy',IyTerm);
    
    obj.domainVariables.set('kappa_loc',kappaTerm);

    % surface variables

    obj.surfaceVariables.model('comp1');
    obj.surfaceVariables.selection.geom('geom', 1);
    obj.surfaceVariables.selection.named('geom_bpeSurface');
    obj.surfaceVariables.label('Variables valid on surface');

    % for Dirichlet BC in potential:
    % surfaceVariables.set('V', 'phi_bpe-projectZetaPotential(phi)', 'Metal - reaction plane potential difference');
    % for Robin BC:
    % surfaceVariables.set('V', '(phi_bpe-phi)*UT', 'Metal - reaction plane potential difference, with dimension');
    obj.surfaceVariables.set('phi_s','phi_bpe');
    obj.surfaceVariables.set('V', '(phi_s-projectReactionPlaneToSurface(phi))*UT', 'Metal - reaction plane potential difference, with dimension');

    % % cathodic currents
    if isempty(cathodicCurrent)
        cathodicCurrent = '0';
    end
    obj.surfaceVariables.set('i_cathodic', cathodicCurrent, 'Cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.');
    obj.surfaceVariables.set('ii_cathodic', cathodicCurrentDimless, 'Dimensionless cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.');


    % % anodic currents    
    if isempty(anodicCurrent)
        anodicCurrent = '0';
    end
    obj.surfaceVariables.set('i_anodic', anodicCurrent, 'Anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.');
    obj.surfaceVariables.set('ii_anodic', anodicCurrentDimless, 'Dimensionless anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.');

    % totalCurrent = strjoin({cathodicCurrent,anodicCurrent},'');
    if isempty(totalCurrent)
        totalCurrent = '0';
    end
    obj.surfaceVariables.set('i_total', totalCurrent, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.');
    obj.surfaceVariables.set('ii_total', totalCurrentDimless, 'Dimensionless total current due to surface reactions. A positive current means charge flux into the electrolyte domain.');

    % QUESTION: should n be includet in exponent???
    for i = 1:obj.nReactions
        obj.surfaceVariables.set(obj.i_id{i}, reactionCurrent{i}, ['current due to surface reaction ', obj.reactions{i}.reaction]);
        obj.surfaceVariables.set(obj.i_dimless_id{i}, reactionCurrentDimless{i}, ['dimensionless current due to surface reaction ', obj.reactions{i}.reaction]);
   end

    % species Flux
    for i = 1:obj.numberOfSpecies
        obj.surfaceVariables.set(obj.N_id{i}, fluxTerm(i));
        obj.surfaceVariables.set(obj.N_dimless_id{i}, fluxTermDimless{i});
    end
    
    obj.bulkVariables.model('comp1');
    obj.bulkVariables.selection.geom('geom', 1);
    obj.bulkVariables.selection.named('geom_bulkBoundary');
    obj.bulkVariables.label('Variables valid on bulk boundary');
    
    obj.upperBoundaryVariables.model('comp1');
    obj.upperBoundaryVariables.selection.geom('geom', 1);
    obj.upperBoundaryVariables.selection.named('geom_upperBoundary');
    obj.upperBoundaryVariables.label('Variables valid on upper boundary');

    obj.electrodeVariables.model('comp1');
    obj.electrodeVariables.selection.geom('geom', 1);
    obj.electrodeVariables.selection.named('geom_electrodes');
    obj.electrodeVariables.label('Variables valid on both working and counter electrode sides');

    obj.weVariables.model('comp1');
    obj.weVariables.selection.geom('geom', 1);
    obj.weVariables.selection.named('geom_workingElectrode');
    obj.weVariables.label('Variables valid  working electrode sides');

    obj.ceVariables.model('comp1');
    obj.ceVariables.selection.geom('geom', 1);
    obj.ceVariables.selection.named('geom_counterElectrode');
    obj.ceVariables.label('Variables valid on counter electrode sides');
    
    obj.weVariables.set('phi_s', 'phi_we+phi_ref');
    obj.ceVariables.set('phi_s', 'phi_ce+phi_ref');
    
    if(obj.explicitElectrodeGeometry)
        obj.electrodeVariables.set('V', '(phi_s-projectReactionPlaneToSurface(phi))*UT', 'Metal - reaction plane potential difference, with dimension');
        % % cathodic currents
        obj.electrodeVariables.set('i_cathodic', cathodicCurrent, 'Cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.');

        % % anodic currents    
        obj.electrodeVariables.set('i_anodic', anodicCurrent, 'Anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.');

        % totalCurrent = strjoin({cathodicCurrent,anodicCurrent},'');
        obj.electrodeVariables.set('i_total', totalCurrent, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.');

        % QUESTION: should n be includet in exponent???
        for i = 1:obj.nReactions
            obj.electrodeVariables.set(obj.i_id{i}, reactionCurrent{i}, ['current due to surface reaction ', obj.reactions{i}.reaction]);
        end

        % species Flux
        for i = 1:obj.numberOfSpecies
            obj.electrodeVariables.set(obj.N_id(i,:), fluxTerm(i));
        end    
    end
    
    
% weak formulation, 2d
    % set discretization
    obj.m.physics('WeakFormulation').prop('ShapeProperty').set('order', '4');
    obj.m.physics('WeakFormulation').prop('ShapeProperty').set('valueType', 'real');
    obj.m.physics('WeakFormulation').prop('ShapeProperty').set('shapeFunctionType', 'shherm');
    
    obj.m.physics('WeakFormulation').field('dimensionless').component([obj.c_id {'phi'}]);
    obj.m.physics('WeakFormulation').feature('wfeq1').set('weak', wEquations);

    
    % Nernst Planck
    for i=1:obj.numberOfSpecies
        % initial values
        obj.m.physics('WeakFormulation').feature('init1').set(obj.c_id{i},obj.c_0_id{i});
        
        % surface flux
        obj.m.physics('WeakFormulation').feature(obj.surfaceFluxBC_id{i}).selection.named('geom_bpeSurface');
        obj.m.physics('WeakFormulation').feature(obj.surfaceFluxBC_id{i}).set('weakExpression',wSurfaceFluxBC{i});
        
        % bulk concentration
        obj.m.physics('WeakFormulation').feature(obj.bulkConcentrationBC_id{i}).selection.named('geom_lateralBoundary');
        obj.m.physics('WeakFormulation').feature(obj.bulkConcentrationBC_id{i}).set('weakExpression',wBulkConcentrationBC{i});
        obj.m.physics('WeakFormulation').feature(obj.bulkConcentrationBC_id{i}).feature([obj.bulkConcentrationBC_id{i},'_aux']).set('fieldVariableName',obj.lambdaCBulk_id{i});
    end
       
    % Poisson
    
    % initial value
    obj.m.physics('WeakFormulation').feature('init1').set('phi', 'phi0');
    
    % electrode potential
    obj.m.physics('WeakFormulation').feature('ElectrodePotentialBC').selection.named('geom_electrodes');
    obj.m.physics('WeakFormulation').feature('ElectrodePotentialBC').set('weakExpression',wElectrodePotentialBC);
    obj.m.physics('WeakFormulation').feature('ElectrodePotentialBC').feature('ElectrodePotentialBC_aux').set('fieldVariableName','lambdaPhiElectrode');
    
    % surface charge density
    obj.m.physics('WeakFormulation').feature('SurfaceChargeDensityBC').selection.named('geom_bpeSurface');
    obj.m.physics('WeakFormulation').feature('SurfaceChargeDensityBC').set('weakExpression',wSurfaceChargeDensityBC);



%     %% Nernst Planck equations
%     fprintf('  Setting up Nernst-Planck system...\n');
% 
%     % settings independent from transformation
%     obj.NernstPlanckEquation.label('Nernst-Planck equation');
%     obj.NernstPlanckEquation.prop('Units').set('DependentVariableQuantity', 'none');
%     obj.NernstPlanckEquation.prop('Units').set('CustomSourceTermUnit', '1'); % really? check this
% 
% %     obj.FluxAtSurface.selection.named('bpeSurface');
%     obj.FluxAtSurface.selection.named('geom_bpeSurface');
%     obj.FluxAtSurface.label('Flux at surface of BPE');
%     
%     obj.FluxAtElectrodes.selection.named('geom_electrodes');
%     obj.FluxAtElectrodes.label('Flux at feeder electrodes');
%     obj.FluxAtElectrodes.active(false);
% 
%     obj.FluxAtBulkBoundary.selection.named('geom_bulkBoundary');
%     obj.FluxAtBulkBoundary.label('Flux from bulk');
% 
%     obj.BulkConcentrationsBC.selection.named('geom_lateralBoundary');
%     obj.BulkConcentrationsBC.label('Fix bulk concentrations');
%     obj.FluxAtBulkBoundary.active(false);
% 
%     fprintf('    Nernst-Planck equation is used in its classical form.\n');
% 
%     % add species concentration variables c
%     obj.NernstPlanckEquation.field('dimensionless').component(obj.c_id);
% 
%     % set up units for variables (concentrations [c] = mol/m^3) and source term
%     % (flux [R] = mol / (m^3 * s) )
%     obj.NernstPlanckEquation.prop('Units').set('CustomDependentVariableUnit', '1');
% 
%     % define governing pde through diffusivities D and preset convection
%     % coefficients determined by species charge and electrical moblity mu
%     obj.NernstPlanckEquation.feature('cfeq1').set('c',pdeNPDiffusionCoefficientC); % needs row array
% 
%     % OUTSOURCE THIS UGLY STUFF TO ANOTHER FILE
% %         pattern = diag(1:numberOfSpecies);
% %         pdeConvectionCoefficientAlphaMatrix = arrayfun( ... 
% %             @(p) iif(  p == 0, @() '0', true, @() pdeNPConvectionCoefficientAlpha{p} ), ...
% %             pattern, 'UniformOutput',false);
% %        temp = strdiag(pdeNPConvectionCoefficientAlpha);
% %         pdeConvectionCoefficientAlphaMatrix = reshape(pdeConvectionCoefficientAlphaMatrix,1,numberOfSpecies^2);
%     pdeNPConvectionCoefficientAlphaXMatrix = flatten(strdiag(pdeNPConvectionCoefficientAlphaX));
%     pdeNPConvectionCoefficientAlphaYMatrix = flatten(strdiag(pdeNPConvectionCoefficientAlphaY));
% 
%     pdeNPConvectionCoefficientAlphaMatrix = [pdeNPConvectionCoefficientAlphaXMatrix,pdeNPConvectionCoefficientAlphaYMatrix];
%     %NernstPlanckEquation.feature('cfeq1').set('al', prepDiagFromString( pdeConvectionCoefficientAlpha, '%i', numberOfSpecies ));
%     obj.NernstPlanckEquation.feature('cfeq1').set('al', pdeNPConvectionCoefficientAlphaMatrix);
%     % 2d modification:
%     % obj.NernstPlanckEquation.feature('cfeq1').set('al', {'z1*phix' 'z1*phiy'; '0' '0'; '0' '0'; 'z2*phix' 'z2*phiy'});
% 
% 
%     % no sources in domain
%     obj.NernstPlanckEquation.feature('cfeq1').set('f', pdeNPSourceTerm);
% 
%     % Nernst-Planck equation boundary conditions
%     % set species flux at bulk to zero
%     % FluxFromBulk = NernstPlanckEquation.create('FluxFromBulk', 'FluxBoundary', 0);
%     % FluxFromBulk.selection.named('bulkBoundary');
%     % FluxFromBulk.label('Flux from bulk');
% 
%     % set species flux at bpe surface due to surface reactions
% %     if( obj.bpeFluxOn )
% %         fprintf('    Surface species flux is activated.\n');
% %         % modifiedFluxTerm = strcat('surfaceFluxRampFactor*(',N_id,')');
%         obj.FluxAtSurface.set('g', pdeNPbpeNeumannBC);
%         % FluxAtBulkBoundary.set('g', cellstr(strcat('-integrateSurface(',pdeNPbpeNeumannBC,')'))); % necessary?
% 
%         % inflow at surface must correspond to outflow at bulk, and vice versa
%         %  FluxAtSurface.set('g', cellstr(N_id));
%         %  FluxAtBulkBoundary.set('g', cellstr(strcat('-integrateSurface(',N_id,')')));
%         % inflow at surface must correspond to outflow at bulk, and vice versa
% 
%         %FluxAtSurface.set('g', zeros(numberOfSpecies,1)); % deactivate surface flux
% %     else
% %        obj.FluxAtSurface.set('g', zeros(obj.numberOfSpecies,1));
% %     end
%     
%     if( obj.explicitElectrodeGeometry )
%         obj.FluxAtElectrodes.set('g', pdeNPbpeNeumannBC);
%         obj.FluxAtElectrodes.active(true);
%     else
%         obj.FluxAtElectrodes.active(false);
%     end
%     
%     % now species flux at bulk
%     obj.FluxAtBulkBoundary.set('g', zeros(obj.numberOfSpecies,1));
% 
%     % fix bulk concentration
%     obj.BulkConcentrationsBC.set('r', pdeNPbulkDirichletBC);
% 
%     % define initial values
%     % modification depending on whether or not Stern layer thickness is part 
%     % of computation domain
%      for i = 1:obj.numberOfSpecies
%         obj.NernstPlanckEquation.feature('init1').set(obj.c_id{i}, pdeNPinitialValues{i});
%      end        
% 
%     %% Poisson equation
%     fprintf('  Setting up Poisson equation...\n');
% 
%     obj.PoissonEquation.label('Poisson equation');
% 
%     % set single variable phi (electric potential)
%     obj.PoissonEquation.field('dimensionless').component({'phi'});
% 
%     % set up units of variables ([phi] = V) and source term ([rho] = C/m^3)
%     obj.PoissonEquation.prop('Units').set('DependentVariableQuantity', 'none');
%     obj.PoissonEquation.prop('Units').set('CustomDependentVariableUnit', '1');
%     obj.PoissonEquation.prop('Units').set('CustomSourceTermUnit', '1');
% 
%     % define governing pde through permeability
%     obj.PoissonEquation.feature('cfeq1').set('c', pdePDiffusionCoefficientC);
% 
%      % prepare charge density term rho = sum rho_i = F sum Z_i c_i 
%     % chargeDensityTerm = [ 'F_const*(', strjoin( cellstr( rcat( ...
%     %     prepIndices('z',numberOfSpecies),...
%     %     repmat('*',numberOfSpecies,1),...
%     %     prepIndices('c',numberOfSpecies) ) )','+' ), ')'];
%     obj.PoissonEquation.feature('cfeq1').set('f', pdePSourceTerm);
% 
%     % Poisson equation boundary conditions
%     
%     % entire surface isolated Neumann BC
%     % n * phi' = 0
%     obj.entireSurfacePotentialGradientBC.selection.named('geom_entireSurface');
%     obj.entireSurfacePotentialGradientBC.set('g', '0');
%     obj.entireSurfacePotentialGradientBC.set('q', '0');
%     obj.entireSurfacePotentialGradientBC.label('Potential gradient at surface, zero gradient = no current, insulator');
%     obj.entireSurfacePotentialGradientBC.active(false);
% 
% %     obj.insulatorSurfacePotentialGradientBC.selection.named('insulatorAdjacentToBpe');
%     obj.insulatorSurfacePotentialGradientBC.selection.named('geom_insulator');
%     obj.insulatorSurfacePotentialGradientBC.set('g', '0');
%     obj.insulatorSurfacePotentialGradientBC.set('q', '0');
%     obj.insulatorSurfacePotentialGradientBC.label('Potential gradient at insulator between bpe and feeder electrodes, zero gradient = no current, insulator');
%     obj.insulatorSurfacePotentialGradientBC.active(false);
% 
%     obj.bulkBoundaryPotentialGradientBC.selection.named('geom_bulkBoundary');
%     obj.bulkBoundaryPotentialGradientBC.set('g', '0');
%     obj.bulkBoundaryPotentialGradientBC.set('q', '0');
%     obj.bulkBoundaryPotentialGradientBC.label('Potential gradient at bulk boundary, zero gradient = no current, insulator');
%     obj.bulkBoundaryPotentialGradientBC.active(false);
% 
%     % set surface charge at boundaries, Robin BC
%     % phi + lambdaS * phi' = phi_bpe, or
%     % C_Stern*(phi_bpe - phi) = eps0*eps* phi' = - sigma
%     % with lambdaS width of Stern layer, sigma surface charge density
%     obj.bpeChargeDensityBC.selection.named('geom_bpeSurface');
%     obj.bpeChargeDensityBC.set('g', pdePbpeRobinBCg);
%     obj.bpeChargeDensityBC.set('q', pdePbpeRobinBCq);
%     obj.bpeChargeDensityBC.label('Charge density at surface');
% 
%     % potential at surface of bpe
%     obj.bpePotentialBC.selection.named('geom_bpeSurface');
%     obj.bpePotentialBC.set('r', pdePbpeDirichletBC);
%     obj.bpePotentialBC.label('potential at surface of BPE');
%     
% 
%     %according initial values:
%     obj.PoissonEquation.feature('init1').set('phi', pdePinitialValues);
% 
%     % chose, which bc to apply to Poisson equation at bpe surface
%     switch obj.bpePoissonEquationBC 
%         case 'Robin'
%             fprintf('    Robin boundary conditions are used...\n');
%             obj.bpeChargeDensityBC.active(true);
%             obj.bpePotentialBC.active(false);
% 
%         case 'Dirichlet'
%             fprintf('    Dirichlet boundary conditions are used...\n');
%             obj.bpeChargeDensityBC.active(false);
%             obj.bpePotentialBC.active(true);
%     end
% 
%     % bulk net charge
%     % bulkChargeDensityBC = PoissonEquation.create('bulkChargeDensityBC', 'FluxBoundary', 0);
%     % bulkChargeDensityBC.selection.named('bulkBoundary');
%     % bulkChargeDensityBC.set('g', '0');
%     % bulkChargeDensityBC.label('Zero net charge density at bulk');
% 
%     % bulk potential
%     obj.bulkPotentialBC.selection.named('geom_bulkBoundary');
%     obj.bulkPotentialBC.set('r', pdePbulkDirichletBC);
%     obj.bulkPotentialBC.label('Bulk potential');
%     obj.bulkPotentialBC.active(false);
%     
%     % electrode potentials
%     obj.electrodePotentialBC.selection.named('geom_electrodes');
%     obj.electrodePotentialBC.set('r', pdePelectrodeDirichletBC);
%     obj.electrodePotentialBC.label('electrode potentials');
%     % initial values, moved into switch fork above
%     % PoissonEquation.feature('init1').set('phi', 'phi_bulk');

    %% Steady-state condition on BPE surface: zero net current
    obj.zeroSurfaceCurrent.label('zeroSurfaceCurrent');
    obj.zeroSurfaceCurrent.name('zeroSurfaceCurrent');
  
%     obj.zeroSurfaceCurrent.field('dimensionless').field('phi_bpe');
%     obj.zeroSurfaceCurrent.field('dimensionless').component({'phi_bpe'});
%     obj.zeroSurfaceCurrent.selection.named('bpeSurface');
    % zeroSurfaceCurrentEquation.selection.named('surfaceNode');

    % switch off standard equation:
%     obj.zeroSurfaceCurrent.feature('dode1').set('f', '0'); % important, otherwise gives wrong results
%     obj.zeroSurfaceCurrent.feature('dode1').set('da', '0');

    % not sure, which formulation is more exact
    obj.zeroSurfaceCurrentEquation.set('name', 'phi_bpe'); 
    obj.zeroSurfaceCurrentEquation.set('equation', 'integrateSurface(i_total)');
    obj.zeroSurfaceCurrentEquation.set('initialValueU', 'phi_bpe_init');
    
    
%     obj.zeroSurfaceCurrent.active(false);

%     obj.zeroSurfaceCurrentEquation.selection.named('bpeSurface');

    % initial value:
%     obj.zeroSurfaceCurrent.feature('init1').set('initialValueU', '0');
end