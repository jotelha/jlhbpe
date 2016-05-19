%%
%   plots are updated to the dataset, whose tag is stored in the string
%   dset
%
function obj = updatePlotsDimensionless(obj,dset)

    fprintf('  Updating 1d plots for dataset "%s"...\n',dset);
    
    obj.m.result('standardPhiPlotGroup').set('data',dset);
    obj.m.result('standardConcentrationsPlotGroup').set('data',dset);
    obj.m.result('logConcentrationsPlotGroup').set('data',dset);
    obj.m.result('multiPurpose1dPlotGroup').set('data',dset);


    %% update dataset for all plot groups
    
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').label('numerical solution');
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('legends', {'numerical'});
%    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').selection.all; % IMPORTANT
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('expr', 'phi');
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('xdata', 'expr'); % get x data from an expression
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('xdataexpr', 'x');
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('xdataunit', '1');
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('legendmethod', 'manual');
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('legend', true);

    % create plot of analytic Poisson-Boltzmann potential distribution
    % phiPBPlotGroup = m.result.duplicate('phiPBPlotGroup','standardPhiPlotGroup');
    % phiPBPlot_id = phiPBPlotGroup.feature.tags; % get the tags of plots under duplicated node
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').label('analytical for z:z electrolyte after PB');
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('legends', {'analytical-for-z:z-electrolyte-after-PB'});
%    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').selection.all; % IMPORTANT
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('expr', 'phi_pb(x,phi_bpe)');
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('xdata', 'expr'); % get x data from an expression
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('xdataexpr', 'x');
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('xdataunit', '1');
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('legendmethod', 'manual');
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('legend', true);
    
    % use arc length instead of expression 'x' for x-axis
    obj.m.result('standardPhiPlotGroup').feature('standardPhiPlot').set('xdata', 'arc'); 
    obj.m.result('standardPhiPlotGroup').feature('phiPBPlot').set('xdata', 'arc'); 
   

    % obj.phiAtSurfacePlotGroup.set('data',dset);
    % obj.logConcentrationsAtSurfacePlotGroup.set('data',dset);

    %% update potential plots
    obj.m.result('standardPhiPlotGroup').label('potential phi');
    obj.m.result('standardPhiPlotGroup').set('ylabel', 'phi / U_T');
    obj.m.result('standardPhiPlotGroup').set('xlabel', 'distance from BPE surface x / L');
    %m.result('standardPhiPlotGroup').set('ylog', true); % logarithmic scale at y axis

    

    %% create concentrations plot

     for i=1:obj.numberOfSpecies
 %       obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).selection.all; %IMPORTANT\
    % use arc length instead of expression 'x' for x-axis
        obj.m.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('xdata', 'arc'); 
        obj.m.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('expr', obj.c_id{i});
%         obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('xdata', 'expr'); % get x data from an expression
%         obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('xdataexpr', 'x');
%         obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('xdataunit', '1');
%         obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('legendmethod', 'manual');
        obj.m.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).label(obj.c_id{i});
%         obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('legends', obj.c_id{i});
%         obj.result('standardConcentrationsPlotGroup').feature(obj.standardConcentrationsPlot_id{i}).set('legend', true);
        
       
     end
    
    obj.m.result('standardConcentrationsPlotGroup').label('species concentrations c');
    obj.m.result('standardConcentrationsPlotGroup').set('ylabel', 'c / c_ref');
    obj.m.result('standardConcentrationsPlotGroup').set('xlabel', 'distance from BPE surface x / L');
    %result('standardConcentrationsPlotGroup').set('ylog', true); % logarithmic scale at y axis

    %% update log plots

     for i=1:obj.numberOfSpecies
%         obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).selection.all; % IMPORTANT
%         obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('expr', obj.c_id{i});
        obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('expr', sprintf('log(%s)',obj.c_id{i}));
        
%         obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('xdata', 'expr'); % get x data from an expression
        % use arc length instead of expression 'x' for x-axis
        obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('xdata', 'arc'); 
        obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('xdataexpr', 'x');
        obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('xdataunit', '1');
%         obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('legendmethod', 'manual');
        obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).label(sprintf('log(%s/c_ref)',obj.c_id{i}));
%         obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('legends',sprintf('log(%s/c_ref)',obj.c_id{i}));
%         obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('legend', true);
        
        % use arc length instead of expression 'x' for x-axis
        obj.m.result('logConcentrationsPlotGroup').feature(obj.logConcentrationsPlot_id{i}).set('xdata', 'arc'); 
     end
    
    obj.m.result('logConcentrationsPlotGroup').label('logarithmic species concentrations c');
    obj.m.result('logConcentrationsPlotGroup').set('ylabel', 'log(c/c_ref)');
    obj.m.result('logConcentrationsPlotGroup').set('xlabel', 'distance from BPE surface x / L');
%     obj.m.result('logConcentrationsPlotGroup').set('ylog', true); % logarithmic scale at y axis

    %% flexible multi purpose plots
    
     for i=1:obj.nMultiPurpose1dPlots
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).set('data', 'parent');
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).set('expr', '0');
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).set('xdata', 'arc'); 
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).set('xdataexpr', 'x');
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).set('xdataunit', '1');
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).label(obj.multiPurpose1dPlot_id{i});
        obj.m.result('multiPurpose1dPlotGroup').feature(obj.multiPurpose1dPlot_id{i}).set('xdata', 'arc'); 
     end
    
    obj.m.result('multiPurpose1dPlotGroup').label('multiPurpose1dPlotGroup');
    obj.m.result('multiPurpose1dPlotGroup').set('ylabel', 'ylabel');
    obj.m.result('multiPurpose1dPlotGroup').set('xlabel', 'x / L');
%     obj.m.result('multiPurpose1dPlotGroup').set('ylog', true); % logarithmic scale at y axis

    %% run plots
    obj.m.result('standardPhiPlotGroup').run();
    obj.m.result('standardConcentrationsPlotGroup').run();
    obj.m.result('logConcentrationsPlotGroup').run();
    obj.m.result('multiPurpose1dPlotGroup').run();


    fprintf('  Finished updating plots.\n');
end