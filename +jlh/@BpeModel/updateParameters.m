function obj = updateParameters(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% prepare replacement patterns
    
    RT_term = 'R_const*T'; 
    UT_term = 'RT/F_const'; % thermal voltage

    fprintf('  Setting global parameters...\n');
    lengthL     = 'lambdaD/epsilon'; % domain depth
%     widthW      = 'lambdaD/gamma'; % domain width
%     widthWbpe   = 'lambdaD/beta'; 
    lambdaS     = 'delta*lambdaD';
   
    % derived expressions

    ionicStrengthSummand = prepTerm('z_id^2*c_bulk_id','z_id','c_bulk_id',obj.z_id,obj.c_bulk_id);
    ionicStrengthTerm = strcat('1/2*(', strjoin( ionicStrengthSummand','+' ), ')');
    debyeLengthTerm = 'sqrt(epsilon0_const*epsilon_r*RT/(2*F_const^2*IonicStrength))';

    wePotentialTerm       = 'phiSymmetryFactor*deltaPhi';
    cePotentialTerm       = '-(1-phiSymmetryFactor)*deltaPhi';

%     wePotentialTerm       = 'phi_bpe+phiSymmetryFactor*deltaPhi';
%     cePotentialTerm       = 'phi_bpe-(1-phiSymmetryFactor)*deltaPhi';
%     
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

    c_ref = 'IonicStrength';
    
    % bulk ionic conductivity
    kappaSummand = prepTerm('z_id^2*c_bulk_id*D_id','z_id','c_bulk_id','D_id',obj.z_id,obj.c_bulk_id,obj.D_id);
    kappaTerm    = ['F_const/RT*(',strcat(strjoin(kappaSummand','+' )),')'];


    %% set gloabl parameters

    % for general pde
    % gammaTerm = '(D%i*C%ix+z%i*D%i*F_const/RT*phix)';
    % sourceTerm = 'C%ix*(D%i*C%ix+z%i*D%i*F_const/RT*phix)';

    for n=1:obj.numberOfSpecies
        obj.m.param.set(obj.D_id{n}, attachUnit(obj.D(n),'m^2/s'), 'Diffusion coefficient');
        obj.m.param.set(obj.z_id{n}, attachUnit(obj.z(n),'1'), 'Charge number');
        obj.m.param.set(obj.c_bulk_id{n}, attachUnit(obj.c_bulk(n),'mol/m^3'), 'Bulk concentration');
    end

% simulation ramping factors
    % all below should be changed to directly insert numeric values from matlab
    obj.m.param.set('surfacePotentialRampFactor', 1);
    obj.m.param.set('surfaceFluxRampFactor', 1);
    obj.m.param.set('chargeDensityRampFactor', 1);
    obj.m.param.set('deltaPhiRampFactor', 1);
    
% dimensionless measures

    obj.m.param.set('epsilon', num2str(obj.epsilon), 'Dimensionless Debye length scale');
    obj.m.param.set('delta', num2str(obj.delta), 'Dimensionless Stern layer width');
%     obj.m.param.set('gamma', num2str(obj.gamma), 'Dimensionless relation betweeb Debye length and domain width');
%     obj.m.param.set('beta', num2str(obj.beta), 'Dimensionless relation betweeb Debye length and BPE width');
    
    obj.m.param.set('h', '1', 'Cell height');

    obj.m.param.set('w', obj.w, 'Cell width');
    obj.m.param.set('w_bpe', obj.w_bpe, 'Width of BPE');
    obj.m.param.set('w_bulkLeft', obj.w_bulkLeft, 'Width of insulating stretch left of system');
    obj.m.param.set('w_bulkRight', obj.w_bulkRight, 'Width of insulating stretch right of system');
  
    obj.m.param.set('w_insulatorLeft', obj.w_insulatorLeft, 'Width of left insulator');
    obj.m.param.set('w_insulatorRight', obj.w_insulatorRight, 'Width of right insulator');
    obj.m.param.set('w_cathode', obj.w_cathode, 'Width of cathode / WE');
    obj.m.param.set('w_anode', obj.w_anode, 'Width of anode / CE');

    obj.m.param.set('w_mesh', obj.w_mesh, 'Mesh chop width');
    obj.m.param.set('nMeshChops', obj.nMeshChops, 'No of mesh chops on BPE');
    obj.m.param.set('nSegmentsPerChop', obj.nSegmentsPerChop, 'No of surface segments per chop');
    obj.m.param.set('extendedDdlFactor', obj.extendedDdlFactor, 'How many debye lengths should mesh refinemend stretch beyond DDL');
    
 
    obj.m.param.set('epsilon_r', num2str(obj.epsilon_r), 'relative permittivity of electrolyte');

    obj.m.param.set('z_ref', max(abs(obj.z)), 'reference charge number for initial value estimate');
    
    obj.m.param.set('kappa_bulk',kappaTerm, 'ionic conductivity in bulk solution'); % dimensional or dimensionless?

    obj.m.param.set('C_Stern', C_Stern, 'Stern layer capacitance, dimensionless'); % Stern layer capacitance

    % m.param.set('phi_bpe',attachUnit(phi_bpe,'V'),'potential at surface of bpe');
    obj.m.param.set('phi_bpe_init',obj.phi_bpe,'dimensionless potential at surface of bpe');
    obj.m.param.set('deltaPhi',deltaPhiTerm,'potential difference between working electrode and counter electrode');
    obj.m.param.set('phiSymmetryFactor',obj.phiSymmetryFactor,'symmetry of applied potential around bpe potential');

    obj.m.param.set('phi_we',wePotentialTerm,'dimensionless potential at working electrode');
    obj.m.param.set('phi_ce',cePotentialTerm,'dimensionless potential at counter electrode');

    obj.m.param.set('phi_bulk','0','bulk electrolyte as reference');
    obj.m.param.set('phi_ref','0','reference potential');

    % dimensional quantities
    obj.m.param.set('IonicStrength', ionicStrengthTerm, 'Ionic Strength');
    obj.m.param.set('c_ref',c_ref, 'reference concentration for dimensionless formulation');
    obj.m.param.set('D_ref',attachUnit(obj.D_ref,'m^2/s'), 'reference diffusivity for dimensionless formulation');
    
    obj.m.param.set('lambdaD', debyeLengthTerm, 'Debye length');
    obj.m.param.set('kappa', 'lambdaD^(-1)', 'reciprocal of Debye length');
    obj.m.param.set('lambdaS', lambdaS); % width of Stern layer, move to parameter settings
    
    obj.m.param.set('L', lengthL, 'Characteristic length, usually multiple of Debye length');
    obj.m.param.set('H', 'h*L', 'Cell depth');

    obj.m.param.set('W', 'w*L', 'Cell width');
    obj.m.param.set('Wbpe', 'w_bpe*L', 'Width of BPE');
    obj.m.param.set('XleftBoundary', '(-w_bpe/2-w_bulkLeft)*L');
    obj.m.param.set('XrightBoundary', '(w_bpe/2+w_bulkRight)*L');
    
    obj.m.param.set('PHI_bpe',['surfacePotentialRampFactor*',attachUnit(obj.PHI_bpe,'V')],'potential at surface of bpe');
    obj.m.param.set('C_Stern_dimensional', C_Stern_dimensional, 'Stern layer capacitance'); % Stern layer capacitance  

    obj.m.param.set('T', attachUnit(obj.T,'K'), 'Temperature'); % move to parameter settings
    obj.m.param.set('RT', RT_term, 'Molar gas constant * Temperature');
    obj.m.param.set('UT', UT_term, 'thermal voltage, R*T/F');
    
    % setup electrochemical surface reaction parameters
    for i=1:obj.nReactions 
%         obj.m.param.set(obj.k0_id{i},attachUnit(obj.reactions{i}.k0,'m/s'),['standard rate constant of ',obj.reactions{i}.reaction]);
%         obj.m.param.set(obj.E0_id{i},attachUnit(obj.reactions{i}.E0,'V'),['reduction potential of ',obj.reactions{i}.reaction]);
        obj.m.param.set(obj.n_id{i},obj.reactions{i}.n,['number of electrons participating in ',obj.reactions{i}.reaction]);
%         obj.m.param.set(obj.beta_id{i},obj.reactions{i}.beta,['symmetry coefficient of ',obj.reactions{i}.reaction]);
    end

    for i=1:obj.numberOfSpecies
        for j = 1:obj.nReactions
            obj.m.param.set(obj.nu_id{i,j},num2str(obj.reactions{j}.flux(i)), sprintf('Flux of species %d due to reaction %s',i, obj.reactions{j}.reaction));
        end
    end
    
    %% custom parameters
    p = fields(obj.customParameters);
    if( ~isempty(p) )
        for i=1:numel(p)
            q = obj.customParameters.(p{i});
            if isstruct(q)
                if isfield(q,'unit') && isnumeric(q.val)
                    val = attachUnit(q.val,q.unit);
                else
                    val = q.val;
                end
                if isfield(q,'desc')
                    obj.m.param.set(p{i},val,q.desc);
                else
                    obj.m.param.set(p{i},val);
                end
            else
                obj.m.param.set(p{i},q);
            end
        end
    end

    
    % smoothen Bpe BC factor
    obj.m.param.set('smootheningFactor', 1);
end