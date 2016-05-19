classdef Constants 
    properties (Constant)
        % Fundamental Constants, from  http://physics.nist.gov/constants
        VacuumPermittivity  = 8.854187817e-12;  %epsilon0, [F/m]
        ElectronCharge      = 1.602176565e-19;  % e, [C]
        FaradayConstant     = 96485.3365;       % F, [C/mol]
        KBoltzmann          = 1.3806488e-23;    % kB, [J/K]
        MolarGasConstant    = 8.3144621;        % R, [J/(mol K)]
        AvogadroConstant    = 6.022140857e23;   % NA, [mol^-1], https://en.wikipedia.org/wiki/Avogadro_constant
        % attributes of water, from http://en.wikipedia.org/wiki/Properties_of_water
        MolarMassOfWater            = 18.01528*1000; % [kg/mol]
        DensityOfWater              = 999.972; % [kg/m^3], rho_H2O
        SpecificHeatCapacityOfWater = 75.375; % [J/(mol*kg)], C_H2O
        ThermalConductivityOfWater  = 0.58; % [W/(m*K)]
%         RelativePermittivityOfWater = 80.18; % at 20 deg C
        RelativePermittivityOfWater = 78.36; % at 25 deg C, Haynes 6-14
        ConductivityOfWater          = 0.055e-2; %[S/m], http://en.wikipedia.org/wiki/Properties_of_water#Electrical_conductivity

        IonicRadiusH2O 				= 0.138e-9; %[m] http://www1.lsbu.ac.uk/water/ionisoh.html
		IonicRadiusOHm 				= 0.11e-9; %[m], http://www1.lsbu.ac.uk/water/ionisoh.html
		MolarMassOHm				= 17.008/1000; %[kg/mol], http://en.wikipedia.org/wiki/Hydroxide
		MolarMassH3Op 				= 19.0232/1000; %[kg/mol], http://en.wikipedia.org/wiki/Hydronium
		IonicRadiusH3Op 			= 0.1e-9; %[m], http://www1.lsbu.ac.uk/water/hydrogen_ions.html
		MolarVolumeH3Op				= 5.4*100^(-3); %[m^3/mol]
        
        DiffusivityH3Op             = 9.31e-9; % [m^2/s]
        DiffusivityOHm              = 1e-0; 
        % source: Krawiec, H.; Vignal, V. & Akid, R.: Numerical modelling of the electrochemical behaviour of 316L stainless steel based upon static and dynamic experimental microcapillary-based techniques. Electrochimica Acta, 53 , 5252-5259 (2008)

       	
		
		% attributes of SDS
		MolarMassSDS 				= 288.372/1000; % [kg/m^3], google
		%IonicRadiusNap 			= 0.116e-9; %[m], http://en.wikipedia.org/wiki/Ionic_radius
		IonicRadiusNap 				= 0.184e-9; %[m], Dey 2012, Aggregation of sodium dodecylsulfate in aqueous nitric acid medium
        MolarMassNa 				= 22.9898/1000; %[kg/mol], google
		MolarMassNap 				= 22.9898/1000; % sane as Na
		IonicRadiusDSm				= 0.402e-9; %[m], same as above
		MolarMassDSm 				= 265.3895/1000 %[kg/mol], http://pubchem.ncbi.nlm.nih.gov/compound/4329331
		
        DiffusivityNap              = 1e-9; % [m^2/s]
        % source: Krawiec, H.; Vignal, V. & Akid, R.: Numerical modelling of the electrochemical behaviour of 316L stainless steel based upon static and dynamic experimental microcapillary-based techniques. Electrochimica Acta, 53 , 5252-5259 (2008)


        % parameters
        T = 298.15; % [K], 25 deg C

        % abbreviations:
        epsilon_0 = jlh.Constants.VacuumPermittivity;
        e = jlh.Constants.ElectronCharge;
        F = jlh.Constants.FaradayConstant;
        kB = jlh.Constants.KBoltzmann;
        R = jlh.Constants.MolarGasConstant;
        f = jlh.Constants.F / (jlh.Constants.R*jlh.Constants.T); % [C/J = 1/V]
        
        epsilon_rw = jlh.Constants.RelativePermittivityOfWater;
    end
end
