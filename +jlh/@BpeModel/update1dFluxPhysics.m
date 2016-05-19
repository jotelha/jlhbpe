function obj = update1dFluxPhysics(obj)
    d = 1; % dimensionality 1
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
%     anodicTerms = prepTerm(...
%         'n_id*F_const*k0_id*(c_R*exp((1-beta_id)*F_const/RT*(V-E0_id)))',...
%         'n_id','k0_id','c_R','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,reductantConcentrations,obj.beta_id,obj.E0_id);
% 29.02.2016, n in exponent
    anodicTerms = prepTerm(...
        'n_id*F_const*k0_id*(c_R*exp(n_id*(1-beta_id)*F_const/RT*(V-E0_id)))',...
        'n_id','k0_id','c_R','beta_id','E0_id',...
        obj.n_id,obj.k0_id,reductantConcentrations,obj.beta_id,obj.E0_id);
    anodicCurrent = strjoin(anodicTerms','+');

    % cathodic currents
%     cathodicTerms = prepTerm(...
%         '-n_id*F_const*k0_id*(c_Ox*exp(-beta_id*F_const/RT*(V-E0_id)))',...
%         'n_id','k0_id','c_Ox','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,oxidantConcentrations,obj.beta_id,obj.E0_id);    
    cathodicTerms = prepTerm(...
        '-n_id*F_const*k0_id*(c_Ox*exp(-n_id*beta_id*F_const/RT*(V-E0_id)))',...
        'n_id','k0_id','c_Ox','beta_id','E0_id',...
        obj.n_id,obj.k0_id,oxidantConcentrations,obj.beta_id,obj.E0_id);    
    cathodicCurrent = strjoin(cathodicTerms','');

    totalCurrent = strjoin({anodicCurrent,cathodicCurrent},'');

    % current per reaction
    % Butler-Volmer, cathodic direction > 0 in form of, as in Bard p. 96
    % i = n*F*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
    % in anodic direction > 0:
%     reactionCurrent = prepTerm(...
%         'n_id*F_const*k0_id*( - c_Ox*exp(-beta_id*F_const/RT*(V-E0_id)) + c_R*exp((1-beta_id)*F_const/RT*(V-E0_id)) ) ',...
%         'n_id','k0_id','c_Ox','c_R','beta_id','E0_id',...
%         obj.n_id,obj.k0_id,oxidantConcentrations,reductantConcentrations,obj.beta_id,obj.E0_id);
     reactionCurrent = prepTerm(...
        'n_id*F_const*k0_id*( c_R*exp(n_id*(1-beta_id)*F_const/RT*(V-E0_id)) - c_Ox*exp(-n_id*beta_id*F_const/RT*(V-E0_id)) ) ',...
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
%     pdeNPConvectionCoefficientAlphaY = prepTerm('z_id*phiy','z_id',obj.z_id);

    pdeNPSourceTerm = strzeros(obj.numberOfSpecies,1); % dropped '

    % for log concentration PDE
    % pdeGammaTerm = prepTerm('-(D_id*C_idx+z_id*D_id*F_const/RT*phix)',...
    %     'D_id','C_id','z_id',D_id,C_id,z_id); % minus convention
    % pdeSourceTerm = prepTerm('C_idx*(D_id*C_idx+z_id*D_id*F_const/RT*phix)',...
    %     'D_id','C_id','z_id',D_id,C_id,z_id);

    % boundary terms NP:
%     pdeNPbulkDirichletBC    = prepTerm('c_bulk_id/c_ref','c_bulk_id',obj.c_bulk_id); % fixed concentration at bulk
    pdeNPbulkDirichletBC    = prepTerm('c_0_id','c_0_id',obj.c_0_id); % fixed concentration at bulk
    pdeNPbpeNeumannBC       = prepTerm('L*N_id/(D_id*c_ref)','N_id','D_id',obj.N_id,obj.D_id); % species flux at surface


    % initial value terms
%     pdeNPinitialValues = prepTerm('c_pb_id(x,phi_bpe)','c_pb_id',obj.c_pb_id);
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
    pdePbulkDirichletBC = '0';
    
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
   
    pdePbpeRobinBCg     = 'C_Stern*phi_s';
    pdePbpeRobinBCq     = 'C_Stern';

    % for delta = 0 => lambdaS = 0, BC reduces to
    %   phi = phi0

    pdePbpeDirichletBC = 'phi_s';

    % for working and counter electrodes
%     pdePelectrodeDirichletBC = 'phi_pb(y,phi_s)';
    pdePelectrodeDirichletBC = 'phi_s';

    % initial values
    pdePinitialValues = 'phi_pb(x,phi_bpe)';
%     pdePinitialValues = 'phi0';
  
    % dimensionless flux N*:
    %   N   = - (D grad c + z F u c grad phi)
    %       = - c_ref D/L (grad* c* + z c* grad* phi*)
    %   N L / (c_ref D_ref) = - D / D_ref (grad* c* + z c* grad* phi*) = N*
    
    nTerm = prepTerm('(-D_id/D_ref*(cx_id+z_id*c_id*phix))','D_id','cx_id','z_id','c_id',obj.D_id, obj.cx_id,obj.z_id,obj.c_id);
%     nyTerm = prepTerm('(-D_id/D_ref*(cy_id+z_id*c_id*phiy))','D_id','cy_id','z_id','c_id',obj.D_id, obj.cy_id,obj.z_id,obj.c_id);

    nTotalTerm = strcat(strjoin(nTerm','+' ));

    iTerm = prepTerm('z_id*nx_id','z_id','nx_id',obj.z_id,obj.nx_id);
%     iySummand = prepTerm('z_id*ny_id','z_id','ny_id',obj.z_id,obj.ny_id);

    iTotalTerm = strcat(strjoin(iTerm','+' ));
%     iyTerm = strcat(strjoin(iySummand','+' ));
    

    % local ionic conductivity
    kappaSummand = prepTerm('z_id^2*c_id*D_id','z_id','c_id','D_id',obj.z_id,obj.c_id,obj.D_id);
    kappaTerm    = ['F_const/RT*(',strcat(strjoin(kappaSummand','+' )),')'];

    c0Term = prepTerm('c_bulk_id/c_ref','c_bulk_id',obj.c_bulk_id);
    phi0Term = '0';
    %% define variables
    
    % domain variables
    
    % concentrations and log concentrations
    obj.domainVariables.model('comp1');
    obj.domainVariables.selection.geom('geom', 1);
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
            obj.domainVariables.set(obj.nx_id{i},nTerm{i});
            obj.domainVariables.set(obj.ix_id{i},iTerm{i});

%             obj.domainVariables.set(obj.ny_id{i},nyTerm{i});
    end
    obj.domainVariables.set('n_total',nTotalTerm);
    obj.domainVariables.set('i_total',iTotalTerm);
%     obj.domainVariables.set('iy',iyTerm);
    
    obj.domainVariables.set('kappa_loc',kappaTerm);

    obj.electrodeVariables.model('comp1');
    obj.electrodeVariables.selection.geom('geom', d-1);
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

    
    obj.FluxAtElectrodes.selection.named('electrodes');
    obj.FluxAtElectrodes.label('Flux at feeder electrodes');

    obj.BulkConcentrationsBC.selection.named('workingElectrode');
    obj.BulkConcentrationsBC.label('Fix inlet concentrations');

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
%     pdeNPConvectionCoefficientAlphaYMatrix = flatten(strdiag(pdeNPConvectionCoefficientAlphaY));

%     pdeNPConvectionCoefficientAlphaMatrix = [pdeNPConvectionCoefficientAlphaXMatrix,pdeNPConvectionCoefficientAlphaYMatrix];
    %NernstPlanckEquation.feature('cfeq1').set('al', prepDiagFromString( pdeConvectionCoefficientAlpha, '%i', numberOfSpecies ));
    obj.NernstPlanckEquation.feature('cfeq1').set('al', pdeNPConvectionCoefficientAlphaXMatrix);
    % 2d modification:
    % obj.NernstPlanckEquation.feature('cfeq1').set('al', {'z1*phix' 'z1*phiy'; '0' '0'; '0' '0'; 'z2*phix' 'z2*phiy'});


    % no sources in domain
    obj.NernstPlanckEquation.feature('cfeq1').set('f', pdeNPSourceTerm);

    % Nernst-Planck equation boundary conditions
    % set species flux at bulk to zero
    % FluxFromBulk = NernstPlanckEquation.create('FluxFromBulk', 'FluxBoundary', 0);
    % FluxFromBulk.selection.named('bulkBoundary');
    % FluxFromBulk.label('Flux from bulk');

    
    % now species flux at bulk
%     obj.FluxAtElectrodes.set('g', zeros(obj.numberOfSpecies,1));
    obj.FluxAtElectrodes.set('g', obj.N_id);


    % fix bulk concentration
    obj.BulkConcentrationsBC.set('r', pdeNPbulkDirichletBC);

    % define initial values
    % modification depending on whether or not Stern layer thickness is part 
    % of computation domain
     for i = 1:obj.numberOfSpecies
        obj.NernstPlanckEquation.feature('init1').set(obj.c_id{i}, pdeNPinitialValues{i});
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


    %according initial values:
    obj.PoissonEquation.feature('init1').set('phi', pdePinitialValues);

    % bulk net charge
    % bulkChargeDensityBC = PoissonEquation.create('bulkChargeDensityBC', 'FluxBoundary', 0);
    % bulkChargeDensityBC.selection.named('bulkBoundary');
    % bulkChargeDensityBC.set('g', '0');
    % bulkChargeDensityBC.label('Zero net charge density at bulk');

    % bulk potential
%     obj.bulkPotentialBC.selection.named('bulkBoundary');
%     obj.bulkPotentialBC.set('r', pdePbulkDirichletBC);
%     obj.bulkPotentialBC.label('Bulk potential');
%     
    % electrode potentials
    obj.electrodePotentialBC.selection.named('electrodes');
    obj.electrodePotentialBC.set('r', pdePelectrodeDirichletBC);
    obj.electrodePotentialBC.label('electrode potentials');
    % initial values, moved into switch fork above
    % PoissonEquation.feature('init1').set('phi', 'phi_bulk');

    %% Steady-state condition on BPE surface: zero net current
    obj.electrodeCurrentEquation.set('name', 'deltaPhi'); 
    obj.electrodeCurrentEquation.set('equation', 'integrateWE(i_total)-i_ref');
    obj.electrodeCurrentEquation.set('initialValueU', '1');

end