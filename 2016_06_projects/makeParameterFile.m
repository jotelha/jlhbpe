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

ionicStrengthSummand = prepTerm('z_id^2*c_bulk_id','z_id','c_bulk_id',m.z_id,m.c_bulk_id);
ionicStrengthTerm = strcat('1/2*(', strjoin( ionicStrengthSummand','+' ), ')');
debyeLengthTerm = 'sqrt(epsilon0_const*epsilon_r*RT/(2*F_const^2*IonicStrength))';

wePotentialTerm       = 'phiSymmetryFactor*DeltaPHI';
cePotentialTerm       = '-(1-phiSymmetryFactor)*DeltaPHI';

%     wePotentialTerm       = 'phi_bpe+phiSymmetryFactor*deltaPhi';
%     cePotentialTerm       = 'phi_bpe-(1-phiSymmetryFactor)*deltaPhi';
%     
deltaPhiTerm          = sprintf('deltaPhiRampFactor*%e',m.deltaPHI);

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
kappaSummand = prepTerm('z_id^2*c_bulk_id*D_id','z_id','c_bulk_id','D_id',m.z_id,m.c_bulk_id,m.D_id);
kappaTerm    = ['F_const/RT*(',strcat(strjoin(kappaSummand','+' )),')'];


%% set gloabl parameters
parameters = containers.Map;

for n=1:m.numberOfSpecies
    parameters(m.D_id{n}) = {attachUnit(m.D(n),'m^2/s'), 'Diffusion coefficient'};
    parameters(m.z_id{n}) = {attachUnit(m.z(n),'1'), 'Charge number'};
    parameters(m.c_bulk_id{n}) ={attachUnit(m.c_bulk(n),'mol/m^3'), 'Bulk concentration'};
end

% simulation ramping factors
% all below should be changed to directly insert numeric values from matlab
parameters('surfacePotentialRampFactor') = 1;
parameters('surfaceFluxRampFactor') = 1;
parameters('chargeDensityRampFactor') = 1;
parameters('deltaPhiRampFactor') = 1;
% smoothen Bpe BC factor
if exist('smootheningFactor','var')
    parameters('smootheningFactor') = smootheningFactor;
else
    parameters('smootheningFactor') = 1;
end

% dimensionless measures

parameters('epsilon') = {num2str(m.epsilon), 'Dimensionless Debye length scale'};
parameters('delta') = {num2str(m.delta), 'Dimensionless Stern layer width'};
% parameters('h') = {'1', 'Cell height'}; % needs to be CHANGED

parameters('w') = {m.w, 'Cell width'};
parameters('w_bpe') = {m.w_bpe, 'Width of BPE'};
parameters('w_bulkLeft') = {m.w_bulkLeft, 'Width of insulating stretch left of system'};
parameters('w_bulkRight') = {m.w_bulkRight, 'Width of insulating stretch right of system'};

% parameters('w_insulatorLeft') = {m.w_insulatorLeft, 'Width of left insulator'};
% parameters('w_insulatorRight') = {m.w_insulatorRight, 'Width of right insulator'};
% parameters('w_cathode') = {m.w_cathode, 'Width of cathode / WE'};
% parameters('w_anode') = {m.w_anode, 'Width of anode / CE'};

% parameters('w_mesh') = {m.w_mesh, 'Mesh chop width'};
% parameters('nMeshChops') = {m.nMeshChops, 'No of mesh chops on BPE'};
% parameters('nSegmentsPerChop') = {m.nSegmentsPerChop, 'No of surface segments per chop'};
% parameters('extendedDdlFactor') = {m.extendedDdlFactor, 'How many debye lengths should mesh refinemend stretch beyond DDL'};

parameters('epsilon_r') = {num2str(m.epsilon_r), 'relative permittivity of electrolyte'};

parameters('z_ref') = {max(abs(m.z)), 'reference charge number for initial value estimate'};

parameters('kappa_bulk') = {kappaTerm, 'ionic conductivity in bulk solution'}; % dimensional or dimensionless?

parameters('C_Stern') = {C_Stern, 'Stern layer capacitance, dimensionless'}; % Stern layer capacitance

% m.param.set('phi_bpe',attachUnit(phi_bpe,'V'),'potential at surface of bpe');
% parameters('phi_bpe_init') = {m.phi_bpe,'dimensionless potential at surface of bpe'};
parameters('deltaPhi') = 'DeltaPHI/UT'; %{deltaPhiTerm,'potential difference between working electrode and counter electrode'};
parameters('phiSymmetryFactor') = {m.phiSymmetryFactor,'symmetry of applied potential around bpe potential'};

parameters('phi_we') = {'PHI_we/UT','dimensionless potential at working electrode'};
parameters('phi_ce') = {'PHI_ce/UT','dimensionless potential at counter electrode'};

parameters('phi_bulk') = {'0','bulk electrolyte as reference'};
parameters('phi_ref') = {'0','reference potential'};

% dimensional quantities
parameters('IonicStrength') =  {ionicStrengthTerm, 'Ionic Strength'};
parameters('c_ref') = { c_ref, 'reference concentration for dimensionless formulation'};
parameters('D_ref') = { attachUnit(m.D_ref,'m^2/s'), 'reference diffusivity for dimensionless formulation'};

parameters('lambdaD') = {debyeLengthTerm, 'Debye length'};
parameters('kappa') = {'lambdaD^(-1)', 'reciprocal of Debye length'};
parameters('lambdaS') = lambdaS; % width of Stern layer, move to parameter settings

parameters('L') = {lengthL, 'Characteristic length, usually multiple of Debye length'};
parameters('H') = {attachUnit(m.H,'m'), 'Cell height'};

parameters('W') = {'w*L', 'Cell width'};
parameters('Wbpe') = {'w_bpe*L', 'Width of BPE'};
parameters('XleftBoundary') = '(-w_bpe/2-w_bulkLeft)*L';
parameters('XrightBoundary') = '(w_bpe/2+w_bulkRight)*L';

parameters('DeltaPHI') = {deltaPhiTerm,'potential difference between working electrode and counter electrode'}; % 'deltaPhi*UT';
parameters('PHI_bpe') = {['surfacePotentialRampFactor*',attachUnit(m.PHI_bpe,'V')],'potential at surface of bpe'};
parameters('PHI_we') = {wePotentialTerm,'dimensional potential at working electrode'}; %'phi_we*UT';
parameters('PHI_ce') = {cePotentialTerm,'dimensional potential at counter electrode'}; %'phi_ce*UT';
parameters('C_Stern_dimensional') = {C_Stern_dimensional, 'Stern layer capacitance'}; % Stern layer capacitance  

parameters('T') = {attachUnit(m.T,'K'), 'Temperature'}; % move to parameter settings
parameters('RT') = {RT_term, 'Molar gas constant * Temperature'};
parameters('UT') = { UT_term, 'thermal voltage, R*T/F'};

% setup electrochemical surface reaction parameters
for i=1:m.nReactions 
%         parameters(m.k0_id{i},attachUnit(m.reactions{i}.k0,'m/s'),['standard rate constant of ',m.reactions{i}.reaction]);
%         parameters(m.E0_id{i},attachUnit(m.reactions{i}.E0,'V'),['reduction potential of ',m.reactions{i}.reaction]);
    parameters(m.n_id{i}) = {m.reactions{i}.n,['number of electrons participating in ',m.reactions{i}.reaction]};
%         parameters(m.beta_id{i},m.reactions{i}.beta,['symmetry coefficient of ',m.reactions{i}.reaction]);
end

for i=1:m.numberOfSpecies
    for j = 1:m.nReactions
        parameters(m.nu_id{i,j}) = {num2str(m.reactions{j}.flux(i)), sprintf('Flux of species %d due to reaction %s',i, m.reactions{j}.reaction)};
    end
end

%% custom parameters
p = fields(m.customParameters);
if( ~isempty(p) )
    for i=1:numel(p)
        q = m.customParameters.(p{i});
        if isstruct(q)
            if isfield(q,'unit') && isnumeric(q.val)
                val = attachUnit(q.val,q.unit);
            else
                val = q.val;
            end
            if isfield(q,'desc')
                parameters(p{i}) = {val,q.desc};
            else
                parameters(p{i}) = val;
            end
        else
            parameters(p{i}) = q;
        end
    end
end

%% save file
parameter = parameters.keys';
value = parameters.values';
hasDesc = cellfun(@(c) iscell(c),value);
comment = cell(numel(parameter),1);
comment(hasDesc) = cellfun(@(c) c{2}, value(hasDesc),'UniformOutput',false);
value(hasDesc) = cellfun(@(c) c{1}, value(hasDesc),'UniformOutput',false);

T = table(parameter,value,comment);
parameterFile = [pwd,'\',m.projectPath,'\parameters.txt'];
files('parameterFile') = parameterFile;
writetable(T,parameterFile,'WriteRowNames',false,'WriteVariableNames',false,'Delimiter',' ');
%% save file records
jlh.hf.saveMapAsTxt(files,'globalFiles.txt');
