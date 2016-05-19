function E_m = calcMixedPotential(obj)
    %  
    % jlh - mixed potential v0.2, Jan 2016
    %

    % note:
    % cathodic currents positive, anodic negative int total current term

    % calculate equilibrium potential by Nernstequation)
%     calcEeq = @(r) r.E0 + obj.RT / (obj.F*r.n) * log( r.cOx / r.cR );
%     Eeq = cellfun(calcEeq,obj.reactions(1:obj.nReactions));
% 
%     calcK0 = @(r) exp(r.n*obj.F/obj.RT*r.E0);
%     K0 = cellfun(calcK0,obj.reactions(1:obj.nReactions));

    i               = cell(obj.nReactions,1);
    i_cathodic      = cell(obj.nReactions,1);
    i_anodic        = cell(obj.nReactions,1);
    i_tot           = @(V) 0;
    i_cathodic_tot  = @(V) 0;
    i_anodic_tot    = @(V) 0;

    % UNCLEAR: how to handle stochiometric coefficients
    % Butler-Volmer, cathodic direction in form of, as in Bard p. 96
    % i = n*F/RT*k0*( cOx * exp( - beta F/RT (E-E0) ) - cR exp( (1-beta) F/RT (E-E0) ) 
    % anodic direction inf form of Newmann 2004, p. 205
    % i = i0 * (exp( alpha_a F eta / RT ) - exp( alpha_c F eta/ RT ))
    
    % jlh 27.2.2015
    %   standard here: 
    %       cathodic current negative
    %       anodic current positive
    for j = 1:obj.nReactions
% terms without n in exponent
%         cathodicTerm    = @(V) obj.reactions{j}.cOx * exp( - obj.reactions{j}.beta * obj.F/obj.RT * (V - obj.reactions{j}.E0) );
%         anodicTerm      = @(V) obj.reactions{j}.cR  * exp( (1-obj.reactions{j}.beta) * obj.F/obj.RT * (V - obj.reactions{j}.E0) );
% terms with n in exponent, 29.02.2016
        cathodicTerm    = @(V) obj.reactions{j}.cOx * exp( - obj.reactions{j}.n*obj.reactions{j}.beta * obj.F/obj.RT * (V - obj.reactions{j}.E0) );
        anodicTerm      = @(V) obj.reactions{j}.cR  * exp( obj.reactions{j}.n*(1-obj.reactions{j}.beta) * obj.F/obj.RT * (V - obj.reactions{j}.E0) );
        
        i_cathodic{j} = @(V) obj.reactions{j}.n * obj.F * obj.reactions{j}.k0 * cathodicTerm(V);
        i_anodic{j} = @(V) obj.reactions{j}.n * obj.F * obj.reactions{j}.k0 * anodicTerm(V);
        i{j} = @(V) obj.reactions{j}.n * obj.F * obj.reactions{j}.k0 * ( ...
            - cathodicTerm(V) + anodicTerm(V) );
        i_cathodic_tot = @(V) i_cathodic_tot(V) + i_cathodic{j}(V);
        i_anodic_tot = @(V) i_anodic_tot(V) + i_anodic{j}(V);
        i_tot = @(V) i_tot(V) + i{j}(V);
    end
    i_cathodic_tot_log = @(V) log(i_cathodic_tot(V));
    i_anodic_tot_log = @(V) log(i_anodic_tot(V));
    % i_tot_log = @(V) i_cathodic_tot_log(V) - i_anodic_tot_log(V);

    % i0 = @(n,ka,kc,cR,cOx,beta) n*F*ka.^beta * kc.^(1-beta) * cR.^beta * cOx.^(1-beta);
    % anodicTerm = @(n,beta,eta) exp((1-beta)*n*F/RT*eta);
    % cathodicTerm = @(n,beta,eta) exp((-beta)*n*F/RT*eta);
    % in = @(eta,cR,cOx,n,beta,ka,kc) i0(n,ka,kc,cR,cOx,beta)*(anodicTerm(n,beta,eta)-cathodicTerm(n,beta,eta));

    %eta = V - E0 = phi_s - phi_l - E_0


    % determine interval
    extractE0 = @(r) r.E0;
    E0 = cellfun(extractE0,obj.reactions(1:obj.nReactions));
    % mixed potential seems not to always lie between the single standard
    % potentials, especially for irreversible reactions:
%     minE0 = min(E0);
%     maxE0 = max(E0);
    minE0 = -2*max(abs(E0));
    maxE0 = 2*max(abs(E0));
    E_m = fzero( i_tot, [minE0,maxE0]);
    % E_m_byLog = fzero( i_tot_log, [minE0,maxE0]);
    
    % use to plot during debugging:
    %   xtest = (minE0:((maxE0-minE0)/100):maxE0);
    %   plot(xtest,i_cathodic_tot_log(xtest),xtest,i_anodic_tot_log(xtest))
    
    % probe for exchange current densities:
    %   obj.reactions{1}.k0*obj.reactions{1}.cR*obj.reactions{1}.n*obj.F
    %   obj.reactions{2}.k0*obj.reactions{2}.cOx*obj.reactions{2}.n*obj.F
    % probe for transfer coefficients
    %   obj.reactions{1}.n*(1-obj.reactions{1}.beta)
    %   obj.reactions{2}.n*obj.reactions{2}.beta
    
    obj.E_m         = E_m;
    obj.i           = i;
    obj.i_anodic    = i_anodic;
    obj.i_cathodic  = i_cathodic;
    obj.i_tot       = i_tot;
    obj.i_anodic_tot = i_anodic_tot;
    obj.i_cathodic_tot = i_cathodic_tot;

    % species flux at mixepd potential
    % r_OHm = 0 - i0OxygenReduction(E_m) / (2*F);
end