function obj = addBpeStudyStepSequence(obj,parameterName,parameterValues,disablePhysics)
    %% attach new study
    fprintf('  Setting up another parametric sweep study...\n');
    
    % parametric ramping
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

    nStudySteps = numel(disablePhysics);
%     disablePhysics = { ...
%         {   'PoissonEquation/bpeChargeDensityBC', ...
%             'NernstPlanckEquation/FluxAtSurface',...
%             'zeroSurfaceCurrent'} ,...
%         {   'NernstPlanckEquation/FluxAtSurface', ...
%             'zeroSurfaceCurrent' }, ...
%         {   'zeroSurfaceCurrent'} ,...
%         {}  };
        
        
    
    studyNum = obj.m.study.size() + 1;
    study_id = sprintf('study%u',studyNum);

    stationaryStudy = obj.m.study.create(study_id);
    
    solNum = obj.m.sol.size() + 1;
    sol_id = sprintf('sol%u',solNum);
    
    obj.m.sol.create(sol_id);
    obj.m.sol(sol_id).study(study_id);
%     sol.attach(study_id);
    obj.m.sol(sol_id).attach(study_id);


    for i=1:nStudySteps
        studyStep_id = sprintf('stationaryStudyStep%d',i);
        compile_id = sprintf('st%d',i);
        variables_id = sprintf('v%d',i);
        solver_id = sprintf('s%d',i);
        store_id = sprintf('su%d',i);
        
        stationaryStudy.create(studyStep_id, 'Stationary');
        
        obj.m.sol(sol_id).create(compile_id, 'StudyStep');
        obj.m.sol(sol_id).create(variables_id, 'Variables');
        obj.m.sol(sol_id).create(solver_id, 'Stationary');
        obj.m.sol(sol_id).create(store_id, 'StoreSolution');
        
%         obj.m.sol(sol_id).feature(solver_id).feature.remove('fcDef');
%         obj.m.sol(sol_id).feature(solver_id).create('fc1', 'FullyCoupled');
%         obj.m.sol(sol_id).feature(solver_id).create('pDef', 'Parametric');
%         obj.m.sol(sol_id).feature(solver_id).create('p1', 'Parametric');

        % what do those stand for:
        obj.m.study(study_id).feature(studyStep_id).set('initstudyhide', 'on');
        obj.m.study(study_id).feature(studyStep_id).set('initsolhide', 'on');
        obj.m.study(study_id).feature(studyStep_id).set('solnumhide', 'on');
        obj.m.study(study_id).feature(studyStep_id).set('notstudyhide', 'on');
        obj.m.study(study_id).feature(studyStep_id).set('notsolhide', 'on');
        obj.m.study(study_id).feature(studyStep_id).set('notsolnumhide', 'on');
        
        obj.m.study(study_id).feature(studyStep_id).set('disabledphysics', disablePhysics{i});
        obj.m.study(study_id).feature(studyStep_id).set('mesh', {'geom' 'standardMesh'});


        obj.m.study(study_id).feature(studyStep_id).set('useparam', useparam);
        obj.m.study(study_id).feature(studyStep_id).set('plistarr', plistarr);
        obj.m.study(study_id).feature(studyStep_id).set('pname', pname);
        obj.m.study(study_id).feature(studyStep_id).set('sweeptype', 'filled');
        obj.m.study(study_id).feature(studyStep_id).set('pcontinuationmode', pcontinuationmode);
        obj.m.study(study_id).feature(studyStep_id).set('useadvanceddisable', true);
        obj.m.study(study_id).feature(studyStep_id).set('showdistribute', true);

        
        obj.m.sol(sol_id).feature(compile_id).set('studystep', studyStep_id);
         
        if( i > 1)
            obj.m.study(study_id).feature(studyStep_id).set('solnum', 'auto');
            obj.m.study(study_id).feature(studyStep_id).set('useinitsol', 'on');
            obj.m.study(study_id).feature(studyStep_id).set('initmethod', 'sol');
            obj.m.study(study_id).feature(studyStep_id).set('initstudy', study_id);
            
            obj.m.sol(sol_id).feature(variables_id).set('initsol', sol_id);
            obj.m.sol(sol_id).feature(variables_id).set('solnum', 'auto');
            obj.m.sol(sol_id).feature(variables_id).set('initmethod', 'sol');
        end
        
     
        obj.m.sol(sol_id).feature(solver_id).set('probesel', 'none');
        obj.m.sol(sol_id).feature(solver_id).set('nonlin', 'on');
        obj.m.sol(sol_id).feature(solver_id).feature('fcDef').set('dtech', 'hnlin');
        obj.m.sol(sol_id).feature(solver_id).feature('fcDef').set('termonres', 'on');
        obj.m.sol(sol_id).feature(solver_id).feature('fcDef').set('niter', '30');
        obj.m.sol(sol_id).feature(solver_id).feature('fcDef').set('ntermauto', 'itertol');
        obj.m.sol(sol_id).feature(solver_id).feature('fcDef').set('minsteph', '1.0E-16');
        obj.m.sol(sol_id).feature(solver_id).feature('pDef').set('sweeptype', 'filled');
        obj.m.sol(sol_id).feature(solver_id).feature('pDef').set('pcontinuationmode', pcontinuationmode);
        obj.m.sol(sol_id).feature(solver_id).feature('pDef').set('plistarr', plistarr);
        obj.m.sol(sol_id).feature(solver_id).feature('pDef').set('pname', pname);
    end 
end