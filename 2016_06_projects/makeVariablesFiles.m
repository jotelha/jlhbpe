import com.comsol.model.*
import com.comsol.model.util.*
import jlh.hf.*
import jlh.*

%% 2d, dimensional

% note concerning comsol normal vectors
% source: https://www.comsol.de/community/forums/general/thread/9119/
% the normal direction are defined wrt the body/domain, it points normally
% "outward" from a closed volume ("up") and if you define a plane you need
% to say where is the geoemtrical entity of higher value (volume for a
% plane, surface for an edge, but a point has a singularity w.r.t normal ;)
% to define the outward ("up") direction ofthe normal.


% prepare replacement patterns
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
fluxDirections = cellfun( @(r) r.flux', m.reactions, 'UniformOutput',false);
fluxDirections = [fluxDirections{:}];
fluxSummandPattern = prepTerm('-(nu_id*i_id/(F_const*n_id))','i_id','n_id',m.i_id,m.n_id);
specificFluxSummand = arrayfun( @(n) jlh.hf.prepTerm(fluxSummandPattern{n},'nu_id',m.nu_id(:,n)),1:m.nReactions,'UniformOutput',false);
specificFluxSummand = [specificFluxSummand{:}];

fluxTerm = cell(m.numberOfSpecies,1);
for i=1:m.numberOfSpecies 
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
fluxTermDimless = prepTerm('L/(D_id*c_ref)*N_id','D_id','N_id',m.D_id,m.N_id);

% IMPORTANT 
% opposite of convention in comsol
%   anodic currents NEGATIVE, cathodic currents POSITIVE

% comsol convention:
%   anodic currents POSITIVE, cathodic currents NEGATIVE

% anodic currents 
anodicTerms = cellfun( @(r) r.anodicCurrentTerm, m.reactions, 'UniformOutput', false);
emptyAnodicTerms = cellfun( @(c) isempty(c), anodicTerms );
anodicCurrent = strjoin(anodicTerms(~emptyAnodicTerms),'+');

% cathodic currents   
cathodicTerms = cellfun( @(r) r.cathodicCurrentTerm, m.reactions, 'UniformOutput', false );
emptyCathodicTerms = cellfun( @(c) isempty(c), cathodicTerms );
cathodicCurrent = strjoin(cathodicTerms(~emptyCathodicTerms),'');

totalCurrent = strjoin({anodicCurrent,cathodicCurrent},'');

% current per reaction
% Butler-Volmer, cathodic direction > 0 in form of, as in Bard p. 96
% i = n*F/RT*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
% in anodic direction > 0:

anodicTerms{emptyAnodicTerms} = '';
cathodicTerms{emptyCathodicTerms} = '';
reactionCurrent = prepTerm('AnodicCathodic','Anodic','Cathodic',anodicTerms',cathodicTerms');
  
% dimensionless currents
reactionCurrentDimless  = prepTerm('L*i_id/(D_ref*c_ref*F_const)','i_id',m.i_id);
anodicCurrentDimless    = 'L*i_anodic/(D_ref*c_ref*F_const)';
cathodicCurrentDimless  = 'L*i_cathodic/(D_ref*c_ref*F_const)';
totalCurrentDimless     = 'L*i_total/(D_ref*c_ref*F_const)';

% dimensionless flux N*:
%   N   = - (D grad c + z F u c grad phi)
%       = - c_ref D/L (grad* c* + z c* grad* phi*)
%   N L / (c_ref D_ref) = - D / D_ref (grad* c* + z c* grad* phi*) = N*

%dimensionless flux: from dimensional flux
% nxTerm = prepTerm('(-cx_id-z_id*c_id*phix)','D_id','cx_id','z_id','c_id',m.D_id, m.cx_id,m.z_id,m.c_id);
% nyTerm = prepTerm('(-cy_id-z_id*c_id*phiy)','D_id','cy_id','z_id','c_id',m.D_id, m.cy_id,m.z_id,m.c_id);
nxTerm = prepTerm('Nx_id*L/(D_id*c_ref)','Nx_id','D_id',m.Nx_id,m.D_id);
nyTerm = prepTerm('Ny_id*L/(D_id*c_ref)','Ny_id','D_id',m.Ny_id,m.D_id);

% dimensional flux: from Tertiary Current Distribution interface
% NxTerm = prepTerm('(D_id*c_ref/L*nx_id)','D_id','nx_id',m.D_id,m.nx_id); % dimensional
% NyTerm = prepTerm('(D_id*c_ref/L*ny_id)','D_id','ny_id',m.D_id,m.ny_id);
NxTerm = prepTerm('domflux.c_idx','c_id',m.c_id); % dimensional
NyTerm = prepTerm('domflux.c_idy','c_id',m.c_id);

ixSummand = prepTerm('z_id*nx_id','z_id','nx_id',m.z_id,m.nx_id);
iySummand = prepTerm('z_id*ny_id','z_id','ny_id',m.z_id,m.ny_id);

ixTerm = strcat(strjoin(ixSummand','+' ));
iyTerm = strcat(strjoin(iySummand','+' ));

IxSummand = prepTerm('F_const*z_id*Ny_id','z_id','Ny_id',m.z_id,m.Nx_id);
IySummand = prepTerm('F_const*z_id*Ny_id','z_id','Ny_id',m.z_id,m.Ny_id);

% IxTerm = 'D_ref*c_ref*F_const/L*ix';
% IyTerm = 'D_ref*c_ref*F_const/L*iy';

IxTerm = strcat(strjoin(IxSummand','+' ));
IyTerm = strcat(strjoin(IySummand','+' ));

% local ionic conductivity
kappaSummand = prepTerm('z_id^2*c_id*D_id','z_id','c_id','D_id',m.z_id,m.c_id,m.D_id);
kappaTerm    = ['F_const/RT*(',strcat(strjoin(kappaSummand','+' )),')'];

c0TermTertiary2d = prepTerm('c_bulk_id','c_bulk_id',m.c_bulk_id);
% c0Term = prepTerm('c_bulk_id/c_ref+smoothenBpeBC(x)*(cInterpolation_id(y)-c_bulk_id/c_ref)','cInterpolation_id','c_bulk_id',m.cInterpolation_id,m.c_bulk_id);

%     phi0Term = '0';
phi0TermTertiary2d = '-(x-(XrightBoundary+XleftBoundary))/W*DeltaPHI';
% phi0Term = 'smoothenBpeBC(x)*phiInterpolation(y)';

% c0TermDilutedSpecies1d = prepTerm('c_bulk_id+smoothenBpeBC(X/L)*(cOHmInterpolation(x)-c_bulk_id)','c_bulk_id',m.c_bulk_id);
c0TermDilutedSpecies1d = prepTerm('c_id_bulk_ddl(X)*exp(-z_id*(phi0-phi_bulk_ddl(X))/UT)','c_id','z_id',m.c_id,m.z_id);

% phi0TermElectrostatics1d = 'smoothenBpeBC(X/L)*phiInterpolation(x)';
phi0TermElectrostatics1d = '(PHI_bpe-phi_bulk_ddl(X))*pbDecayFunction(x)+phi_bulk_ddl(X)';

%% Diluted Species and Electrostatics, 1d

domainVariables = containers.Map;
surfaceVariables = containers.Map;
bulkVariables = containers.Map;
weVariables = containers.Map;
ceVariables = containers.Map;

domainVariables('phi0') = phi0TermElectrostatics1d;
    
for i = 1:m.numberOfSpecies
   domainVariables(m.c_0_id{i}) = c0TermDilutedSpecies1d{i};
   domainVariables(m.nx_id{i}) = nxTerm{i};
   domainVariables(m.Nx_id{i}) = NxTerm{i};
end

domainVariables('ix') = ixTerm;
domainVariables('Ix') = IxTerm;
domainVariables('kappa_loc') = kappaTerm;

% surface variable
surfaceVariables('phi_s') = 'PHI_bpe';
surfaceVariables('V') = { 'phi_s-phi_bulk_ddl(X)', 'Metal - reaction plane potential difference, with dimension'};

surfaceVariables('i_cathodic') = {cathodicCurrent, 'Cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.'};
surfaceVariables('ii_cathodic') = {cathodicCurrentDimless, 'Dimensionless cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.'};

% % anodic currents    
surfaceVariables('i_anodic') = {anodicCurrent, 'Anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.'};
surfaceVariables('ii_anodic') = {anodicCurrentDimless, 'Dimensionless anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.'};

% totalCurrent = strjoin({cathodicCurrent,anodicCurrent},'');
surfaceVariables('i_total') = {totalCurrent, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.'};
surfaceVariables('ii_total') = {totalCurrentDimless, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.'};

for i = 1:m.nReactions
    surfaceVariables(m.i_id{i}) = {reactionCurrent{i}, ['current due to surface reaction ', m.reactions{i}.reaction]};
    surfaceVariables(m.i_dimless_id{i}) = {reactionCurrentDimless{i}, ['dimensionless current due to surface reaction ', m.reactions{i}.reaction]};
end

% species Flux
for i = 1:m.numberOfSpecies
    surfaceVariables(m.N_id{i}) = fluxTerm{i};
    surfaceVariables(m.N_dimless_id{i}) = fluxTermDimless{i};
end

bulkVariables('phi_s') = 'phi_bulk';

files('domainVariablesDilutedSpeciesAndElectrostatics1d') = [pwd,'\',m.projectPath,'\domainVariablesDilutedSpeciesAndElectrostatics1d.txt'];
files('surfaceVariablesDilutedSpeciesAndElectrostatics1d') = [pwd,'\',m.projectPath,'\surfaceVariablesDilutedSpeciesAndElectrostatics1d.txt'];
files('bulkVariablesDilutedSpeciesAndElectrostatics1d') = [pwd,'\',m.projectPath,'\bulkVariablesDilutedSpeciesAndElectrostatics1d.txt'];
writeParameterTextFile(domainVariables,files('domainVariablesDilutedSpeciesAndElectrostatics1d'));
writeParameterTextFile(surfaceVariables,files('surfaceVariablesDilutedSpeciesAndElectrostatics1d'));
writeParameterTextFile(bulkVariables,files('bulkVariablesDilutedSpeciesAndElectrostatics1d'));
    
%% Tertiary Current Distribution, 2d

% domain variables

% domainVariables('domainVariables').model('comp1');
% domainVariables('domainVariables').selection.geom('geom', 2);
% domainVariables('domainVariables').selection.all;

domainVariables('phi0') = phi0TermTertiary2d;
for i = 1:m.numberOfSpecies
    domainVariables(m.c_0_id{i}) = c0TermTertiary2d{i};
end

% domainVariables('domainVariables').label('Variables valid over entire domain');

for i = 1:m.numberOfSpecies
        domainVariables(m.nx_id{i}) = nxTerm{i};
        domainVariables(m.ny_id{i}) = nyTerm{i};
        domainVariables(m.Nx_id{i}) = NxTerm{i};
        domainVariables(m.Ny_id{i}) = NyTerm{i};
end
domainVariables('ix') = ixTerm;
domainVariables('iy') = iyTerm;
domainVariables('Ix') = IxTerm;
domainVariables('Iy') = IyTerm;

domainVariables('kappa_loc') = kappaTerm;

% surface variables

% for Dirichlet BC in potential:
% m.variable('surfaceVariables').set('V', 'phi_bpe-projectZetaPotential(phi)', 'Metal - reaction plane potential difference');
% for Robin BC:
% m.variable('surfaceVariables').set('V', '(phi_bpe-phi)*UT', 'Metal - reaction plane potential difference, with dimension');
surfaceVariables('phi_s') = 'PHI_bpe';
surfaceVariables('V') = {'phi_s-phi', 'Metal - reaction plane potential difference, with dimension'};

% % cathodic currents
if isempty(cathodicCurrent)
    cathodicCurrent = '0';
end
surfaceVariables('i_cathodic') = {cathodicCurrent, 'Cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.'};
surfaceVariables('ii_cathodic') = {cathodicCurrentDimless, 'Dimensionless cathodic current due to surface reactions. A positive cathodic current means charge flux into the electrolyte domain.'};

% % anodic currents    
if isempty(anodicCurrent)
    anodicCurrent = '0';
end
surfaceVariables('i_anodic') = {anodicCurrent, 'Anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.'};
surfaceVariables('ii_anodic') = {anodicCurrentDimless, 'Dimensionless anodic current due to surface reactions. A positive anodic current means charge flux out of the electrolyte domain.'};

% totalCurrent = strjoin({cathodicCurrent,anodicCurrent},'');
if isempty(totalCurrent)
    totalCurrent = '0';
end
surfaceVariables('i_total') = {totalCurrent, 'Total current due to surface reactions. A positive current means charge flux into the electrolyte domain.'};
surfaceVariables('ii_total') = {totalCurrentDimless, 'Dimensionless total current due to surface reactions. A positive current means charge flux into the electrolyte domain.'};

for i = 1:m.nReactions
    surfaceVariables(m.i_id{i}) = {reactionCurrent{i}, ['current due to surface reaction ', m.reactions{i}.reaction]};
    surfaceVariables(m.i_dimless_id{i}) = {reactionCurrentDimless{i}, ['dimensionless current due to surface reaction ', m.reactions{i}.reaction]};
end

% species Flux
for i = 1:m.numberOfSpecies
    surfaceVariables(m.N_id{i}) = fluxTerm{i};
    surfaceVariables(m.N_dimless_id{i}) = fluxTermDimless{i};
end

% electrode variables

weVariables('phi_s') = 'PHI_we';
ceVariables('phi_s') = 'PHI_ce';

files('domainVariablesTertiaryCurrentDistribution2d') = [pwd,'\',m.projectPath,'\domainVariablesTertiaryCurrentDistribution2d.txt'];
files('surfaceVariablesTertiaryCurrentDistribution2d') = [pwd,'\',m.projectPath,'\surfaceVariablesTertiaryCurrentDistribution2d.txt'];
files('weVariablesTertiaryCurrentDistribution2d') = [pwd,'\',m.projectPath,'\weVariablesTertiaryCurrentDistribution2d.txt'];
files('ceVariablesTertiaryCurrentDistribution2d') = [pwd,'\',m.projectPath,'\ceVariablesTertiaryCurrentDistribution2d.txt'];
writeParameterTextFile(domainVariables,files('domainVariablesTertiaryCurrentDistribution2d'));
writeParameterTextFile(surfaceVariables,files('surfaceVariablesTertiaryCurrentDistribution2d'));
writeParameterTextFile(weVariables,files('weVariablesTertiaryCurrentDistribution2d'));
writeParameterTextFile(ceVariables,files('ceVariablesTertiaryCurrentDistribution2d'));

%% Diluted Species and Electrostatics, 2d
surfaceVariables('V') = {'phi_s-extrudeZetaPlaneToSurface(phi)', 'Metal - reaction plane potential difference, with dimension'};

files('domainVariablesDilutedSpeciesAndElectrostatics2d') = [pwd,'\',m.projectPath,'\domainVariablesDilutedSpeciesAndElectrostatics2d.txt'];
files('surfaceVariablesDilutedSpeciesAndElectrostatics2d') = [pwd,'\',m.projectPath,'\surfaceVariablesDilutedSpeciesAndElectrostatics2d.txt'];
files('weVariablesDilutedSpeciesAndElectrostatics2d') = [pwd,'\',m.projectPath,'\weVariablesDilutedSpeciesAndElectrostatics2d.txt'];
files('ceVariablesDilutedSpeciesAndElectrostatics2d') = [pwd,'\',m.projectPath,'\ceVariablesDilutedSpeciesAndElectrostatics2d.txt'];
writeParameterTextFile(domainVariables,files('domainVariablesDilutedSpeciesAndElectrostatics2d'));
writeParameterTextFile(surfaceVariables,files('surfaceVariablesDilutedSpeciesAndElectrostatics2d'));
writeParameterTextFile(weVariables,files('weVariablesDilutedSpeciesAndElectrostatics2d'));
writeParameterTextFile(ceVariables,files('ceVariablesDilutedSpeciesAndElectrostatics2d'));

%% save file records
jlh.hf.saveMapAsTxt(files,'globalFiles.txt');