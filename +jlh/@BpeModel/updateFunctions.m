function obj = updateFunctions(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% create replacement patterns
    % analytic Poisson-Boltzmann solution, uses concentration of first species
    % assuming symmetric z:z electrolyte
    % phiPBFunction = '4*RT/(abs(z1)*F_const)*atanh( tanh( abs(z1)*F_const * phi_bpe/(4*RT) )*exp(-x/lambdaD))';
    % phiPBxFunction = '- 4*RT/(lambdaD*abs(z1)*F_const)* sinh( abs(z1) * F_const* phi_bpe/(4*RT) )* cosh( abs(z1) * F_const* phi_bpe/(4*RT) ) / ( exp(x/lambdaD)* cosh( abs(z1) * F_const* phi_bpe/(4*RT) )^2 - exp(-x/lambdaD) * sinh( abs(z1) * F_const* phi_bpe/(4*RT) )^2)';
    % correction for Stern layer:
    % phiPBFunction = '4*RT/(abs(z1)*F_const)*atanh( tanh( abs(z1)*F_const * phi_bpe/(4*RT) )*exp(-(x+lambdaS)/lambdaD))';
    % dimensionless: 
    %   phiPB* = phiPB/UT = F/RT*phiPB
    %   phiPBx* = L/UT * phiPBx
    phiPBFunction = '4/abs(z_ref)*atanh( tanh( abs(z_ref)*phi0/4 )*exp(-(x/epsilon+delta)))'; % dimensionless
    % phiPBxFunction = '- 4*RT/(lambdaD*abs(z_ref)*F_const)* sinh( abs(z_ref) * F_const* phi_bpe/(4*RT) )* cosh( abs(z_ref) * F_const* phi_bpe/(4*RT) ) / ( exp((x+lambdaS)/lambdaD)* cosh( abs(z_ref) * F_const* phi_bpe/(4*RT) )^2 - exp(-(x+lambdaS)/lambdaD) * sinh( abs(z_ref) * F_const* phi_bpe/(4*RT) )^2)';
    phiPBxFunction = '- 4/(epsilon*abs(z_ref))* sinh( abs(z_ref) * phi0/4 )* cosh( abs(z_ref) * phi0/4 ) / ( exp(x/epsilon+delta)* cosh( abs(z_ref) * phi0/4 )^2 - exp(-(x/epsilon+delta)) * sinh( abs(z_ref) * phi0/4 )^2)';


    % cPBFunction = prepTerm('c_bulk_id*exp(-z_id*F_const*phi_pb(x)/RT)',...
    %     'c_bulk_id','z_id',c_bulk_id,z_id);
    % dimensionless:
    %   cPB* = c_bulk/c_ref * exp( -z*phiPB*)
    cPBFunction = prepTerm('c_bulk_id/c_ref*exp(-z_id*phi_pb(x,phi0))',...
        'c_bulk_id','z_id',obj.c_bulk_id,obj.z_id);
    
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
        obj.c_pb{i}.set('funcname', obj.c_pb_id{i});
        obj.c_pb{i}.set('args', {'x' 'phi0'});
        obj.c_pb{i}.set('expr', cPBFunction{i});
        obj.c_pb{i}.set('fununit', '1');
        obj.c_pb{i}.set('argunit', '1');
        obj.c_pb{i}.set('plotargs', {'x' '0' '1'; 'phi0' 'phi_bpe' 'phi_bpe'});
        %c_pb{i}.label('Concentration distribution in symmetric binary electrolyte due to Poisson-Boltzmann model');
    end
    
    % smoothening function
    obj.m.func('smoothenBpeBC').set('upper', 'w_bpe/2');
    obj.m.func('smoothenBpeBC').set('smooth', 'epsilon*smootheningFactor');
    obj.m.func('smoothenBpeBC').set('funcname', 'smoothenBpeBC');
    obj.m.func('smoothenBpeBC').set('lower', '-w_bpe/2');

end