function obj = addParametricSweepStudy(obj,parameterName,parameterValues)
    %% attach new study
    fprintf('  Setting up another parametric sweep study...\n');

    studyNum = obj.m.study.size() + 1;
    study_id = sprintf('study%u',studyNum);

    obj.stationaryStudy = obj.m.study.create(study_id);
    obj.stationaryStudyStep1 = obj.stationaryStudy.create('stationaryStudyStep1', 'Stationary');

    % what do those stand for:
    obj.stationaryStudyStep1.set('initstudyhide', 'on');
    obj.stationaryStudyStep1.set('initsolhide', 'on');
    obj.stationaryStudyStep1.set('solnumhide', 'on');
    obj.stationaryStudyStep1.set('notstudyhide', 'on');
    obj.stationaryStudyStep1.set('notsolhide', 'on');
    obj.stationaryStudyStep1.set('notsolnumhide', 'on');
    obj.stationaryStudyStep1.set('mesh', {'geom' 'standardMesh'});


    solNum = obj.m.sol.size() + 1;
    sol_id = sprintf('sol%u',solNum);

    sol = obj.m.sol.create(sol_id);
    sol.study(study_id);
    sol.attach(study_id);
    sol.create('st1', 'StudyStep');
    sol.create('v1', 'Variables');
    obj.Stationary = sol.create('s1', 'Stationary');
    obj.FullyCoupled = obj.Stationary.create('fc1', 'FullyCoupled');

    obj.Stationary.set('nonlin', 'on'); % for nonlinear problems, I am not sure about its indications


    %% solver settings
    obj.FullyCoupled.set('dtech','hnlin'); % treat as highly non-linear problem
    %obj.FullyCoupled.set('dtech', 'const'); % constant damping factor
    %obj.FullyCoupled.set('initsteph','1e-4'); % initial damping factor
    obj.FullyCoupled.set('minsteph','1e-16'); % minimal damping factor
    %obj.FullyCoupled.set('ntol','1e-4'); % tolerance 
    %obj.FullyCoupled.set('ntolfact','1'); % tolerance factor, multiplied with main solver tolerance
    %obj.FullyCoupled.set('useratelimit','off'); % disable stopping due to slow convergence
    obj.FullyCoupled.set('termonres','on'); % 'on' terminates on relative residual, off terminates solution-based

    obj.FullyCoupled.set('ntermauto', 'itertol'); % termination criterion: iteration or tolerance
    obj.FullyCoupled.set('niter', '30'); % terminate afte 5iterations

    %% parametric ramping
    fprintf('      Surface potential is ramped.\n');
    pname = parameterName;
    plistarr = cellstr( cellfun(@(c) num2str(c), parameterValues,'UniformOutput',false));
    punit = '';
    useparam = 'on';
    pcontinuationmode = 'no';
    % OR
%     pcontinuationmode = 'manual'; % or 'auto', 'last'
%     pcontinuation = pname;
%     pout = 'plist'; % or 'psteps'

    obj.stationaryStudyStep1.set('pname', pname);
    obj.stationaryStudyStep1.set('plistarr', plistarr);
    obj.stationaryStudyStep1.set('punit', punit);
    obj.stationaryStudyStep1.set('useparam', useparam);
    obj.stationaryStudyStep1.set('pcontinuationmode', pcontinuationmode);
%     obj.stationaryStudyStep1.set('activate',...
%         {'NernstPlanckEquation' 'on' 'PoissonEquation' 'on' 'zeroSurfaceCurrent' 'off'});
    obj.stationaryStudyStep1.set('sweeptype', 'filled');
    %obj.stationaryStudyStep1.set('plot', 'on');

    obj.Parametric = obj.Stationary.create('p1', 'Parametric');
    obj.Stationary.feature.remove('fcDef'); % what is this for?
    obj.Parametric.set('punit', punit);
    obj.Parametric.set('plistarr', plistarr);
    obj.Parametric.set('pname', pname);
    obj.Parametric.set('pcontinuationmode', pcontinuationmode); % deactivate continuation
    obj.Parametric.set('sweeptype', 'filled'); % all combinations
end