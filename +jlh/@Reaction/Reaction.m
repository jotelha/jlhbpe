classdef Reaction < handle
    properties (SetAccess = public)
        name
        reaction
        oxidants % names of oxidants
        oxidantConcentrations % variable names of oxidant concentrations
        cOx % numeric value of bulk Oxidant concentration term
        nuOxidants % stochiometric coefficients 
        reductants
        reductantConcentrations
        cR
        nuReductants 
        n % number of electrons participating
        beta % symmetry coefficient of reaction
        k0 % standard rate constant of reaction
        E0 % standard reduction potential
        flux % stochiometric coefficients in array
    end
end