function obj = makeChoppedGeometry(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% modifying geometry
    fprintf('  Modifying geometry...\n');

%     obj.geom.repairTol(1e-16); % should be smaller than smallest elements
%     % obj.space.set('intervals','many');
%     % space.set('p','0,lambdaD,L'); % an interval with three points
%     % space.set('p','0,lambdaD/L,1'); % an interval with three points. dimensionless
% 
%     %electrode_surface = space.
% 
%     obj.ddl.set('size', {'w_mesh' 'epsilon'});
%     obj.ddl.set('pos', {'-w_bpe/2' '0'});
% %     obj.space.setIndex('layer', 'epsilon', 0);
% %     obj.space.set('layername', {'dl'}); % diffuse layer
%     obj.space.set('size', {'w_mesh' '1-epsilon'});
%     obj.space.set('pos', {'-w_bpe/2' 'epsilon'});  
%     
%     obj.spaceArray.selection('input').set({'space','ddl'});
%     obj.spaceArray.set('displ', 'w_mesh 0');
%     obj.spaceArray.set('size', 'nMeshChops 1');
%     
%     obj.gapLeft.set('size', {'w_bulkLeft' '1'});
%     obj.gapLeft.set('pos', {'-w_bpe/2-w_bulkLeft' '0'});
%     
%     obj.gapRight.set('size', {'w_bulkRight' '1'});
%     obj.gapRight.set('pos', {'w_bpe/2' '0'});
%     
%     obj.leftEdgeOfBpe.setIndex('p', '-w_bpe/2', 0, 0);
%     obj.leftEdgeOfBpe.setIndex('p', '0', 1, 0);
%     obj.rightEdgeOfBpe.setIndex('p', 'w_bpe/2', 0, 0);
%     obj.rightEdgeOfBpe.setIndex('p', '0', 1, 0);
%     
%     if(obj.explicitElectrodeGeometry)
%         obj.leftEdgeOfCathode.setIndex('p', '-w_bpe/2-w_insulatorLeft-w_cathode', 0, 0);
%         obj.leftEdgeOfCathode.setIndex('p', '0', 1, 0);
%         obj.rightEdgeOfCathode.setIndex('p', '-w_bpe/2-w_insulatorLeft', 0, 0);
%         obj.rightEdgeOfCathode.setIndex('p', '0', 1, 0);
%         
%         obj.leftEdgeOfAnode.setIndex('p', 'w_bpe/2+w_insulatorRight', 0, 0);
%         obj.leftEdgeOfAnode.setIndex('p', '0', 1, 0);
%         obj.rightEdgeOfAnode.setIndex('p', 'w_bpe/2+w_insulatorRight+w_anode', 0, 0);
%         obj.rightEdgeOfAnode.setIndex('p', '0', 1, 0);
%     end
%     
%     obj.geom.feature('fin').set('repairtol', '1.0E-10');
%     % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
%     %obj.geom.feature('fin').set('repairtol', '1.0E-16');
% 	obj.geom.run('fin');
% %     obj.geom.run;
obj.geom = obj.m.geom.create('geom',2); % 1 for one dimension

obj.m.geom('geom').repairTol(1.0E-10);

%% cumulative selections
obj.m.geom('geom').selection.create('upperBoundaryCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('upperBoundaryCumulative').label('upperBoundaryCumulative');
obj.m.geom('geom').selection.create('leftBoundarySelection', 'CumulativeSelection');
obj.m.geom('geom').selection('leftBoundarySelection').label('leftBoundarySelection');
obj.m.geom('geom').selection.create('meshChopPrototype', 'CumulativeSelection');
obj.m.geom('geom').selection('meshChopPrototype').label('meshChopPrototype');
obj.m.geom('geom').selection.create('rightBoundarySelection', 'CumulativeSelection');
obj.m.geom('geom').selection('rightBoundarySelection').label('rightBoundarySelection');
obj.m.geom('geom').selection.create('gapLeftSelection', 'CumulativeSelection');
obj.m.geom('geom').selection('gapLeftSelection').label('gapLeftSelection');
obj.m.geom('geom').selection.create('gapRightSelection', 'CumulativeSelection');
obj.m.geom('geom').selection('gapRightSelection').label('gapRightSelection');
obj.m.geom('geom').selection.create('zetaPlaneCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('zetaPlaneCumulative').label('zetaPlaneCumulative');
obj.m.geom('geom').selection.create('bpeSurfaceCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('bpeSurfaceCumulative').label('bpeSurfaceCumulative');
obj.m.geom('geom').selection.create('entireSurfaceCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('entireSurfaceCumulative').label('entireSurfaceCumulative');
obj.m.geom('geom').selection.create('regionOfFirstDebyeLength', 'CumulativeSelection');
obj.m.geom('geom').selection('regionOfFirstDebyeLength').label('regionOfFirstDebyeLength');

obj.m.geom('geom').selection.create('meshChopArrayCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('meshChopArrayCumulative').label('meshChopArrayCumulative');
obj.m.geom('geom').selection.create('domainUnionCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('domainUnionCumulative').label('domainUnionCumulative');

obj.m.geom('geom').selection.create('unitCellPrototypeCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('unitCellPrototypeCumulative').label('unitCellPrototypeCumulative');
obj.m.geom('geom').selection.create('unitCellPrototypeBoundaryCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('unitCellPrototypeBoundaryCumulative').label('unitCellPrototypeBoundaryCumulative');

obj.m.geom('geom').selection.create('eastwardThinningCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('eastwardThinningCumulative').label('eastwardThinningCumulative');
obj.m.geom('geom').selection.create('westwardThinningCumulative', 'CumulativeSelection');
obj.m.geom('geom').selection('westwardThinningCumulative').label('westwardThinningCumulative');



%% geometry
% requires at least 4 bpe mesh subsections

% space prototype for replication
obj.m.geom('geom').create('space', 'Rectangle');
obj.m.geom('geom').feature('space').set('contributeto', 'meshChopPrototype');
obj.m.geom('geom').feature('space').set('createselection', 'on');
obj.m.geom('geom').feature('space').set('size', {'w_mesh' '1-epsilon'});
obj.m.geom('geom').feature('space').set('pos', {'-w_bpe/2' 'epsilon'});
% obj.m.geom('geom').feature('space').set('pos', {'-w_bpe/2+w_mesh' 'epsilon'});

% selections on space prototype
obj.m.geom('geom').create('zetaPlaneEdgePrototype', 'ExplicitSelection');
obj.m.geom('geom').feature('zetaPlaneEdgePrototype').set('contributeto', 'zetaPlaneCumulative');
obj.m.geom('geom').feature('zetaPlaneEdgePrototype').label('zetaPlaneEdgePrototype');
obj.m.geom('geom').feature('zetaPlaneEdgePrototype').selection('selection').init(1);
obj.m.geom('geom').feature('zetaPlaneEdgePrototype').selection('selection').set('space(1)', [1]);
obj.m.geom('geom').create('upperBoundaryEdgePrototype', 'ExplicitSelection');
obj.m.geom('geom').feature('upperBoundaryEdgePrototype').set('contributeto', 'upperBoundaryCumulative');
obj.m.geom('geom').feature('upperBoundaryEdgePrototype').label('upperBoundaryEdgePrototype');
obj.m.geom('geom').feature('upperBoundaryEdgePrototype').selection('selection').init(1);
obj.m.geom('geom').feature('upperBoundaryEdgePrototype').selection('selection').set('space(1)', [3]);

% ddl prototype for replication
obj.m.geom('geom').create('ddl', 'Rectangle');
obj.m.geom('geom').feature('ddl').set('contributeto', 'meshChopPrototype');
obj.m.geom('geom').feature('ddl').set('createselection', 'on');
obj.m.geom('geom').feature('ddl').set('size', {'w_mesh' 'epsilon'});
obj.m.geom('geom').feature('ddl').set('pos', {'-w_bpe/2' '0'});
% obj.m.geom('geom').feature('ddl').set('pos', {'-w_bpe/2+w_mesh' '0'});
obj.m.geom('geom').create('bpeSurfaceEdgePrototype', 'ExplicitSelection');
obj.m.geom('geom').feature('bpeSurfaceEdgePrototype').set('contributeto', 'bpeSurfaceCumulative');
obj.m.geom('geom').feature('bpeSurfaceEdgePrototype').label('bpeSurfaceEdgePrototype');
obj.m.geom('geom').feature('bpeSurfaceEdgePrototype').selection('selection').init(1);
obj.m.geom('geom').feature('bpeSurfaceEdgePrototype').selection('selection').set('ddl(1)', [1]);
obj.m.geom('geom').create('entrieSurfaceEdgePrototype', 'ExplicitSelection');
obj.m.geom('geom').feature('entrieSurfaceEdgePrototype').set('contributeto', 'entireSurfaceCumulative');
obj.m.geom('geom').feature('entrieSurfaceEdgePrototype').label('entrieSurfaceEdgePrototype');
obj.m.geom('geom').feature('entrieSurfaceEdgePrototype').selection('selection').init(1);
obj.m.geom('geom').feature('entrieSurfaceEdgePrototype').selection('selection').set('ddl(1)', [1]);
obj.m.geom('geom').create('regionOfFirstDebyeLengthPrototype', 'ExplicitSelection');
obj.m.geom('geom').feature('regionOfFirstDebyeLengthPrototype').set('contributeto', 'regionOfFirstDebyeLength');
obj.m.geom('geom').feature('regionOfFirstDebyeLengthPrototype').selection('selection').set('ddl(1)', [1]);

obj.m.geom('geom').create('extendedDdl', 'Rectangle');
obj.m.geom('geom').feature('extendedDdl').set('contributeto', 'meshChopPrototype');
obj.m.geom('geom').feature('extendedDdl').set('createselection', 'on');
obj.m.geom('geom').feature('extendedDdl').set('size', {'w_mesh' 'extendedDdlFactor*epsilon'});
obj.m.geom('geom').feature('extendedDdl').set('pos', {'-w_bpe/2' '0'});

% create the ends of domain next to bpe
obj.m.geom('geom').create('gapLeft', 'Rectangle');
obj.m.geom('geom').feature('gapLeft').set('contributeto', 'gapLeftSelection');
obj.m.geom('geom').feature('gapLeft').set('createselection', 'on');
obj.m.geom('geom').feature('gapLeft').set('size', {'w_bulkLeft' '1'});
obj.m.geom('geom').feature('gapLeft').set('pos', {'-w_bpe/2-w_bulkLeft' '0'});

obj.m.geom('geom').create('leftGapEntireSurfaceEdgeContribution', 'ExplicitSelection');
obj.m.geom('geom').feature('leftGapEntireSurfaceEdgeContribution').set('contributeto', 'entireSurfaceCumulative');
obj.m.geom('geom').feature('leftGapEntireSurfaceEdgeContribution').label('leftGapEntireSurfaceEdgeContribution');
obj.m.geom('geom').feature('leftGapEntireSurfaceEdgeContribution').selection('selection').init(1);
obj.m.geom('geom').feature('leftGapEntireSurfaceEdgeContribution').selection('selection').set('gapLeft(1)', [1]);
obj.m.geom('geom').create('leftBoundaryEdge', 'ExplicitSelection');
obj.m.geom('geom').feature('leftBoundaryEdge').set('contributeto', 'leftBoundarySelection');
obj.m.geom('geom').feature('leftBoundaryEdge').label('leftBoundaryEdge');
obj.m.geom('geom').feature('leftBoundaryEdge').selection('selection').init(1);
obj.m.geom('geom').feature('leftBoundaryEdge').selection('selection').set('gapLeft(1)', [4]);
obj.m.geom('geom').create('leftGapUpperBoundaryEdgeContribution', 'ExplicitSelection');
obj.m.geom('geom').feature('leftGapUpperBoundaryEdgeContribution').set('contributeto', 'upperBoundaryCumulative');
obj.m.geom('geom').feature('leftGapUpperBoundaryEdgeContribution').label('leftGapUpperBoundaryEdgeContribution');
obj.m.geom('geom').feature('leftGapUpperBoundaryEdgeContribution').selection('selection').init(1);
obj.m.geom('geom').feature('leftGapUpperBoundaryEdgeContribution').selection('selection').set('gapLeft(1)', [3]);

obj.m.geom('geom').create('leftBoundaryOfSurface', 'ExplicitSelection');
obj.m.geom('geom').feature('leftBoundaryOfSurface').label('leftBoundaryOfSurface');
obj.m.geom('geom').feature('leftBoundaryOfSurface').selection('selection').init(0);
obj.m.geom('geom').feature('leftBoundaryOfSurface').selection('selection').set('gapLeft(1)', [1]);
obj.m.geom('geom').create('leftEdgeOfBpe', 'ExplicitSelection');
obj.m.geom('geom').feature('leftEdgeOfBpe').label('leftEdgeOfBpe');
obj.m.geom('geom').feature('leftEdgeOfBpe').selection('selection').init(0);
obj.m.geom('geom').feature('leftEdgeOfBpe').selection('selection').set('gapLeft(1)', [2]);

% right end piece
obj.m.geom('geom').create('gapRight', 'Rectangle');
obj.m.geom('geom').feature('gapRight').set('contributeto', 'gapRightSelection');
obj.m.geom('geom').feature('gapRight').set('createselection', 'on');
obj.m.geom('geom').feature('gapRight').set('size', {'w_bulkRight' '1'});
obj.m.geom('geom').feature('gapRight').set('pos', {'w_bpe/2' '0'});

obj.m.geom('geom').create('rightGapEntireSurfaceEdgeContribution', 'ExplicitSelection');
obj.m.geom('geom').feature('rightGapEntireSurfaceEdgeContribution').set('contributeto', 'entireSurfaceCumulative');
obj.m.geom('geom').feature('rightGapEntireSurfaceEdgeContribution').label('rightGapEntireSurfaceEdgeContribution');
obj.m.geom('geom').feature('rightGapEntireSurfaceEdgeContribution').selection('selection').init(1);
obj.m.geom('geom').feature('rightGapEntireSurfaceEdgeContribution').selection('selection').set('gapRight(1)', [1]);
obj.m.geom('geom').create('rightBoundaryEdge', 'ExplicitSelection');
obj.m.geom('geom').feature('rightBoundaryEdge').set('contributeto', 'rightBoundarySelection');
obj.m.geom('geom').feature('rightBoundaryEdge').label('rightBoundaryEdge');
obj.m.geom('geom').feature('rightBoundaryEdge').selection('selection').init(1);
obj.m.geom('geom').feature('rightBoundaryEdge').selection('selection').set('gapRight(1)', [2]);
obj.m.geom('geom').create('rightGapUpperBoundaryEdgeContribution', 'ExplicitSelection');
obj.m.geom('geom').feature('rightGapUpperBoundaryEdgeContribution').set('contributeto', 'upperBoundaryCumulative');
obj.m.geom('geom').feature('rightGapUpperBoundaryEdgeContribution').label('rightGapUpperBoundaryEdgeContribution');
obj.m.geom('geom').feature('rightGapUpperBoundaryEdgeContribution').selection('selection').init(1);
obj.m.geom('geom').feature('rightGapUpperBoundaryEdgeContribution').selection('selection').set('gapRight(1)', [3]);

obj.m.geom('geom').create('rightBoundaryOfSurface', 'ExplicitSelection');
obj.m.geom('geom').feature('rightBoundaryOfSurface').label('rightBoundaryOfSurface');
obj.m.geom('geom').feature('rightBoundaryOfSurface').selection('selection').init(0);
obj.m.geom('geom').feature('rightBoundaryOfSurface').selection('selection').set('gapRight(1)', [2]);
obj.m.geom('geom').create('rightEdgeOfBpe', 'ExplicitSelection');
obj.m.geom('geom').feature('rightEdgeOfBpe').label('rightEdgeOfBpe');
obj.m.geom('geom').feature('rightEdgeOfBpe').selection('selection').init(0);
obj.m.geom('geom').feature('rightEdgeOfBpe').selection('selection').set('gapRight(1)', [1]);

% replicate bpe pieces
obj.m.geom('geom').create('spaceArray', 'Array');
obj.m.geom('geom').feature('spaceArray').set('size', 'nMeshChops 1');
obj.m.geom('geom').feature('spaceArray').set('displ', {'w_mesh' '0'});
obj.m.geom('geom').feature('spaceArray').selection('input').set({'space' 'ddl' 'extendedDdl'});
obj.m.geom('geom').feature('spaceArray').set('contributeto', 'meshChopArrayCumulative');

% obj.m.geom('geom').create('leftEdgeOfBpe', 'Point');
% obj.m.geom('geom').feature('leftEdgeOfBpe').active(false);
% obj.m.geom('geom').feature('leftEdgeOfBpe').setIndex('p', '-w_bpe/2', 0, 0);
% obj.m.geom('geom').feature('leftEdgeOfBpe').setIndex('p', '0', 1, 0);
% obj.m.geom('geom').create('rightEdgeOfBpe', 'Point');
% obj.m.geom('geom').feature('rightEdgeOfBpe').setIndex('p', 'w_bpe/2', 0, 0);
% obj.m.geom('geom').feature('rightEdgeOfBpe').setIndex('p', '0', 1, 0);
% obj.m.geom('geom').create('leftEdgeOfCathode', 'Point');
% obj.m.geom('geom').feature('leftEdgeOfCathode').active(false);
% obj.m.geom('geom').create('rightEdgeOfCathode', 'Point');
% obj.m.geom('geom').feature('rightEdgeOfCathode').active(false);
% obj.m.geom('geom').create('leftEdgeOfAnode', 'Point');
% obj.m.geom('geom').feature('leftEdgeOfAnode').active(false);
% obj.m.geom('geom').create('rightEdgeOfAnode', 'Point');
% obj.m.geom('geom').feature('rightEdgeOfAnode').active(false);

% rectangle for mesh refinement at edges of bpe
obj.m.geom('geom').create('eastwardThinningInner', 'Rectangle');
obj.m.geom('geom').feature('eastwardThinningInner').label('eastwardThinningInner');
obj.m.geom('geom').feature('eastwardThinningInner').set('contributeto', 'eastwardThinningCumulative');
obj.m.geom('geom').feature('eastwardThinningInner').set('createselection', 'on');
obj.m.geom('geom').feature('eastwardThinningInner').set('size', {'epsilon' 'epsilon'});
obj.m.geom('geom').feature('eastwardThinningInner').set('pos', {'w_bpe/2-epsilon' '0'});

obj.m.geom('geom').create('westwardThinningInner', 'Rectangle');
obj.m.geom('geom').feature('westwardThinningInner').label('westwardThinningInner');
obj.m.geom('geom').feature('westwardThinningInner').set('contributeto', 'westwardThinningCumulative');
obj.m.geom('geom').feature('westwardThinningInner').set('createselection', 'on');
obj.m.geom('geom').feature('westwardThinningInner').set('size', {'epsilon' 'epsilon'});
obj.m.geom('geom').feature('westwardThinningInner').set('pos', {'-w_bpe/2' '0'});

obj.m.geom('geom').create('eastwardThinningOuter', 'Rectangle');
obj.m.geom('geom').feature('eastwardThinningOuter').label('eastwardThinningOuter');
obj.m.geom('geom').feature('eastwardThinningOuter').set('contributeto', 'eastwardThinningCumulative');
obj.m.geom('geom').feature('eastwardThinningOuter').set('createselection', 'on');
obj.m.geom('geom').feature('eastwardThinningOuter').set('size', {'epsilon' 'epsilon'});
obj.m.geom('geom').feature('eastwardThinningOuter').set('pos', {'-w_bpe/2-epsilon' '0'});

obj.m.geom('geom').create('westwardThinningOuter', 'Rectangle');
obj.m.geom('geom').feature('westwardThinningOuter').label('westwardThinningOuter');
obj.m.geom('geom').feature('westwardThinningOuter').set('contributeto', 'westwardThinningCumulative');
obj.m.geom('geom').feature('westwardThinningOuter').set('createselection', 'on');
obj.m.geom('geom').feature('westwardThinningOuter').set('size', {'epsilon' 'epsilon'});
obj.m.geom('geom').feature('westwardThinningOuter').set('pos', {'w_bpe/2' '0'});

% % union of all domains
% obj.m.geom('geom').create('domainUnion', 'UnionSelection');
% obj.m.geom('geom').feature('domainUnion').set('entitydim', '-1'); % object
% obj.m.geom('geom').feature('domainUnion').label('domainUnion');
% obj.m.geom('geom').feature('domainUnion').set('input', {'meshChopArrayCumulative',...
%     'gapLeftSelection','gapRightSelection',...
%     'westwardThinningCumulative','eastwardThinningCumulative'});
% % obj.m.geom('geom').feature('domainUnion').set('contributeto', 'domainUnionCumulative');
% 
% 
% obj.m.geom('geom').create('uniteDomains','Union');
% % % obj.m.geom('geom').feature('uniteDomains').set('entitydim', '2');
% obj.m.geom('geom').feature('uniteDomains').selection('input').named('domainUnion');
% % %     'gapLeftSelection','gapRightSelection',...
% % %     'westwardThinningCumulative','eastwardThinningCumulative'});
% % obj.m.geom('geom').feature('uniteDomains').('meshChopArray');

% unions and assemblies
obj.m.geom('geom').create('allSelection', 'BoxSelection');
obj.m.geom('geom').feature('allSelection').set('entitydim', '-1');
obj.m.geom('geom').feature('allSelection').label('allSelection');

obj.m.geom('geom').create('ddlSelection', 'UnionSelection');
obj.m.geom('geom').feature('ddlSelection').set('entitydim', '-1');
obj.m.geom('geom').feature('ddlSelection').label('dllSelection');
obj.m.geom('geom').feature('ddlSelection').set('input', {'westwardThinningOuter' 'eastwardThinningOuter' ...
    'westwardThinningInner' 'eastwardThinningInner' 'ddl'});

obj.m.geom('geom').create('ddlUnion', 'Union');
obj.m.geom('geom').feature('ddlUnion').label('ddlUnion');
obj.m.geom('geom').feature('ddlUnion').set('selresult', 'on');
obj.m.geom('geom').feature('ddlUnion').selection('input').named('ddlSelection');

obj.m.geom('geom').create('allExceptDdl', 'Difference');
obj.m.geom('geom').feature('allExceptDdl').label('allExceptDdl');
obj.m.geom('geom').feature('allExceptDdl').set('keep', true);
obj.m.geom('geom').feature('allExceptDdl').set('selresult', 'on');
obj.m.geom('geom').feature('allExceptDdl').selection('input').named('allSelection');
obj.m.geom('geom').feature('allExceptDdl').selection('input2').named('ddlUnion');

obj.m.geom('geom').create('toDelete', 'DifferenceSelection');
obj.m.geom('geom').feature('toDelete').set('entitydim', '-1');
obj.m.geom('geom').feature('toDelete').label('toDelete');
obj.m.geom('geom').feature('toDelete').set('add', {'allSelection'});
obj.m.geom('geom').feature('toDelete').set('subtract', {'allExceptDdl' 'ddlUnion'});

obj.m.geom('geom').create('delete', 'Delete');
obj.m.geom('geom').feature('delete').selection('input').init;
obj.m.geom('geom').feature('delete').selection('input').named('toDelete');

% points
obj.m.geom('geom').create('leftBoundaryOfZetaPlanePoint', 'Point');
obj.m.geom('geom').feature('leftBoundaryOfZetaPlanePoint').label('leftBoundaryOfZetaPlanePoint');
obj.m.geom('geom').feature('leftBoundaryOfZetaPlanePoint').setIndex('p', '-w_bpe/2-w_bulkLeft', 0, 0);
obj.m.geom('geom').feature('leftBoundaryOfZetaPlanePoint').setIndex('p', 'epsilon', 1, 0);
obj.leftBoundaryOfZetaPlane = obj.m.geom('geom').create('leftBoundaryOfZetaPlane', 'ExplicitSelection');
obj.m.geom('geom').feature('leftBoundaryOfZetaPlane').label('leftBoundaryOfZetaPlane');
obj.m.geom('geom').feature('leftBoundaryOfZetaPlane').selection('selection').init(0);
obj.m.geom('geom').feature('leftBoundaryOfZetaPlane').selection('selection').set('leftBoundaryOfZetaPlanePoint(1)', [1]);
obj.m.geom('geom').create('rightBoundaryOfZetaPlanePoint', 'Point');
obj.m.geom('geom').feature('rightBoundaryOfZetaPlanePoint').label('rightBoundaryOfZetaPlanePoint');
obj.m.geom('geom').feature('rightBoundaryOfZetaPlanePoint').setIndex('p', 'w_bpe/2+w_bulkRight', 0, 0);
obj.m.geom('geom').feature('rightBoundaryOfZetaPlanePoint').setIndex('p', 'epsilon', 1, 0);
obj.rightBoundaryOfZetaPlane = obj.m.geom('geom').create('rightBoundaryOfZetaPlane', 'ExplicitSelection');
obj.m.geom('geom').feature('rightBoundaryOfZetaPlane').label('rightBoundaryOfZetaPlane');
obj.m.geom('geom').feature('rightBoundaryOfZetaPlane').selection('selection').init(0);
obj.m.geom('geom').feature('rightBoundaryOfZetaPlane').selection('selection').set('rightBoundaryOfZetaPlanePoint(1)', [1]);
obj.insulatorAdjacentToBpe = obj.m.geom('geom').create('insulatorAdjacentToBpe', 'DifferenceSelection');
obj.m.geom('geom').feature('insulatorAdjacentToBpe').set('entitydim', '1');
obj.m.geom('geom').feature('insulatorAdjacentToBpe').label('insulatorAdjacentToBpe');
obj.m.geom('geom').feature('insulatorAdjacentToBpe').set('add', {'entireSurfaceCumulative'});
obj.m.geom('geom').feature('insulatorAdjacentToBpe').set('subtract', {'bpeSurfaceCumulative'});
obj.lateralBoundary = obj.m.geom('geom').create('lateralBoundary', 'UnionSelection');
obj.m.geom('geom').feature('lateralBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('lateralBoundary').label('lateralBoundary');
obj.m.geom('geom').feature('lateralBoundary').set('input', {'leftBoundarySelection' 'rightBoundarySelection'});
obj.counterElectrode = obj.m.geom('geom').create('workingElectrode', 'UnionSelection');
obj.m.geom('geom').feature('workingElectrode').set('entitydim', '1');
obj.m.geom('geom').feature('workingElectrode').label('workingElectrode 1');
obj.m.geom('geom').feature('workingElectrode').set('input', {'leftBoundarySelection'});
obj.counterElectrode = obj.m.geom('geom').create('counterElectrode', 'UnionSelection');
obj.m.geom('geom').feature('counterElectrode').set('entitydim', '1');
obj.m.geom('geom').feature('counterElectrode').label('counterElectrode');
obj.m.geom('geom').feature('counterElectrode').set('input', {'rightBoundarySelection'});
obj.electrodes = obj.m.geom('geom').create('electrodes', 'UnionSelection');
obj.m.geom('geom').feature('electrodes').set('entitydim', '1');
obj.m.geom('geom').feature('electrodes').label('electrodes');
obj.m.geom('geom').feature('electrodes').set('input', {'rightBoundarySelection' 'leftBoundarySelection'});
obj.bulkBoundary = obj.m.geom('geom').create('bulkBoundary', 'UnionSelection');
obj.m.geom('geom').feature('bulkBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('bulkBoundary').label('bulkBoundary');
obj.m.geom('geom').feature('bulkBoundary').set('input', {'rightBoundarySelection' 'leftBoundarySelection' 'upperBoundaryCumulative'});

obj.bpeSurface = obj.m.geom('geom').create('bpeSurface', 'UnionSelection');
obj.m.geom('geom').feature('bpeSurface').set('entitydim', '1');
obj.m.geom('geom').feature('bpeSurface').label('bpeSurface');
obj.m.geom('geom').feature('bpeSurface').set('input', {'bpeSurfaceCumulative'});
obj.zetaPlane = obj.m.geom('geom').create('zetaPlane', 'UnionSelection');
obj.m.geom('geom').feature('zetaPlane').set('entitydim', '1');
obj.m.geom('geom').feature('zetaPlane').label('zetaPlane');
obj.m.geom('geom').feature('zetaPlane').set('input', {'zetaPlaneCumulative'});
obj.entireSurface = obj.m.geom('geom').create('entireSurface', 'UnionSelection');
obj.m.geom('geom').feature('entireSurface').set('entitydim', '1');
obj.m.geom('geom').feature('entireSurface').label('entireSurface');
obj.m.geom('geom').feature('entireSurface').set('input', {'entireSurfaceCumulative'});
obj.upperBoundary = obj.m.geom('geom').create('upperBoundary', 'UnionSelection');
obj.m.geom('geom').feature('upperBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('upperBoundary').label('upperBoundary');
obj.m.geom('geom').feature('upperBoundary').set('input', {'upperBoundaryCumulative'});
obj.insulator = obj.m.geom('geom').create('insulator', 'UnionSelection');
obj.m.geom('geom').feature('insulator').set('entitydim', '1');
obj.m.geom('geom').feature('insulator').label('insulator');
obj.m.geom('geom').feature('insulator').set('input', {'insulatorAdjacentToBpe'});

%% selections important for unit cell meshing

% take second piece
obj.m.geom('geom').create('unitCellPrototype', 'BoxSelection');
obj.m.geom('geom').feature('unitCellPrototype').label('unitCellPrototype');
obj.m.geom('geom').feature('unitCellPrototype').set('xmin', '-w_bpe/2+3*w_mesh/2');
obj.m.geom('geom').feature('unitCellPrototype').set('xmax', '-w_bpe/2+3*w_mesh/2');
obj.m.geom('geom').feature('unitCellPrototype').set('ymin', '-Inf');
obj.m.geom('geom').feature('unitCellPrototype').set('ymax', 'Inf');
obj.m.geom('geom').feature('unitCellPrototype').set('condition', 'intersects');
obj.m.geom('geom').feature('unitCellPrototype').set('contributeto', 'unitCellPrototypeCumulative');

obj.m.geom('geom').create('unitCellPrototypeUpperDomain', 'BoxSelection');
obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').label('unitCellPrototypeUpperDomain');
obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('xmin', '-w_bpe/2+3*w_mesh/2');
obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('xmax', '-w_bpe/2+3*w_mesh/2');
% obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('ymin', 'epsilon/2+1/2');
obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('ymin', '3*epsilon/2'); % to include extended BC
obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('ymax', 'epsilon/2+1/2');
obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('condition', 'intersects');
% obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('contributeto', 'unitCellPrototypeCumulative');

obj.m.geom('geom').create('unitCellPrototypeLowerDomain', 'BoxSelection');
obj.m.geom('geom').feature('unitCellPrototypeLowerDomain').label('unitCellPrototypeLowerDomain');
obj.m.geom('geom').feature('unitCellPrototypeLowerDomain').set('xmin', '-w_bpe/2+3*w_mesh/2');
obj.m.geom('geom').feature('unitCellPrototypeLowerDomain').set('xmax', '-w_bpe/2+3*w_mesh/2');
obj.m.geom('geom').feature('unitCellPrototypeLowerDomain').set('ymin', 'epsilon/2');
obj.m.geom('geom').feature('unitCellPrototypeLowerDomain').set('ymax', 'epsilon/2');
obj.m.geom('geom').feature('unitCellPrototypeLowerDomain').set('condition', 'intersects');
% obj.m.geom('geom').feature('unitCellPrototypeUpperDomain').set('contributeto', 'unitCellPrototypeCumulative');


% obj.m.geom('geom').create('unitCellPrototypeBoundary', 'BoxSelection');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').label('unitCellPrototypeBoundary');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('entitydim', '1');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('xmax', '-w_bpe/2+w_mesh');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('xmin', '-w_bpe/2');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('ymin', '0');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('ymax', '1');
% obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('contributeto', 'unitCellPrototypeBoundaryCumulative');
obj.m.geom('geom').create('unitCellPrototypeBoundary', 'AdjacentSelection');
obj.m.geom('geom').feature('unitCellPrototypeBoundary').label('unitCellPrototypeBoundary');
obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('entitydim', '2');
obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('input', {'unitCellPrototypeCumulative'});
obj.m.geom('geom').feature('unitCellPrototypeBoundary').set('contributeto', 'unitCellPrototypeBoundaryCumulative');

obj.m.geom('geom').create('unitCellPrototypeBpeSurface', 'IntersectionSelection');
obj.m.geom('geom').feature('unitCellPrototypeBpeSurface').label('unitCellPrototypeBpeSurface');
obj.m.geom('geom').feature('unitCellPrototypeBpeSurface').set('entitydim', '1');
obj.m.geom('geom').feature('unitCellPrototypeBpeSurface').set('input', {'unitCellPrototypeBoundaryCumulative', 'bpeSurfaceCumulative'});

obj.m.geom('geom').create('unitCellPrototypeZetaPlane', 'IntersectionSelection');
obj.m.geom('geom').feature('unitCellPrototypeZetaPlane').label('unitCellPrototypeZetaPlane');
obj.m.geom('geom').feature('unitCellPrototypeZetaPlane').set('entitydim', '1');
obj.m.geom('geom').feature('unitCellPrototypeZetaPlane').set('input', {'unitCellPrototypeBoundaryCumulative', 'zetaPlaneCumulative'});

obj.m.geom('geom').create('unitCellPrototypeUpperBoundary', 'IntersectionSelection');
obj.m.geom('geom').feature('unitCellPrototypeUpperBoundary').label('unitCellPrototypeUpperBoundary');
obj.m.geom('geom').feature('unitCellPrototypeUpperBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('unitCellPrototypeUpperBoundary').set('input', {'unitCellPrototypeBoundaryCumulative', 'upperBoundaryCumulative'});

obj.m.geom('geom').create('boundariesAdjacentToBpeSurface', 'AdjacentSelection');
obj.m.geom('geom').feature('boundariesAdjacentToBpeSurface').label('boundariesAdjacentToBpeSurface');
obj.m.geom('geom').feature('boundariesAdjacentToBpeSurface').set('entitydim', '1');
obj.m.geom('geom').feature('boundariesAdjacentToBpeSurface').set('input', {'bpeSurface'});

obj.m.geom('geom').create('boundariesAdjacentToUpperBoundary', 'AdjacentSelection');
obj.m.geom('geom').feature('boundariesAdjacentToUpperBoundary').label('boundariesAdjacentToUpperBoundary');
obj.m.geom('geom').feature('boundariesAdjacentToUpperBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('boundariesAdjacentToUpperBoundary').set('input', {'upperBoundary'});

obj.m.geom('geom').create('unitCellPrototypeLowerLateralBoundary', 'IntersectionSelection');
obj.m.geom('geom').feature('unitCellPrototypeLowerLateralBoundary').label('unitCellPrototypeLowerLateralBoundary');
obj.m.geom('geom').feature('unitCellPrototypeLowerLateralBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('unitCellPrototypeLowerLateralBoundary').set('input', {'unitCellPrototypeBoundaryCumulative', 'boundariesAdjacentToBpeSurface'});

obj.m.geom('geom').create('unitCellPrototypeUpperLateralBoundary', 'IntersectionSelection');
obj.m.geom('geom').feature('unitCellPrototypeUpperLateralBoundary').label('unitCellPrototypeUpperLateralBoundary');
obj.m.geom('geom').feature('unitCellPrototypeUpperLateralBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('unitCellPrototypeUpperLateralBoundary').set('input', {'unitCellPrototypeBoundaryCumulative', 'boundariesAdjacentToUpperBoundary'});

% obj.m.geom('geom').create('boundariesAdjacentToEntireSurface', 'AdjacentSelection');
% obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').label('boundariesAdjacentToEntireSurface');
% obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('entitydim', '1');
% obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('input', {'entireSurface'});

obj.m.geom('geom').create('boundariesAdjacentToEntireSurface', 'BoxSelection');
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').label('boundariesAdjacentToEntireSurface');
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('entitydim',1);
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('xmin', '-Inf');
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('xmax', 'Inf');
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('ymin', 'epsilon/2');
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('ymax', 'epsilon/2');
obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('condition', 'intersects');
% obj.m.geom('geom').feature('boundariesAdjacentToEntireSurface').set('contributeto', 'boundariesAdjacentToEntireSurfaceCumulative');

% edge distribution selections
obj.m.geom('geom').create('westwardThinningEdges', 'DifferenceSelection');
obj.m.geom('geom').feature('westwardThinningEdges').label('westwardThinningEdges');
obj.m.geom('geom').feature('westwardThinningEdges').set('entitydim', '1');
obj.m.geom('geom').feature('westwardThinningEdges').set('add', {'westwardThinningCumulative'});
obj.m.geom('geom').feature('westwardThinningEdges').set('subtract', {'boundariesAdjacentToEntireSurface'});

obj.m.geom('geom').create('eastwardThinningEdges', 'DifferenceSelection');
obj.m.geom('geom').feature('eastwardThinningEdges').label('eastwardThinningEdges');
obj.m.geom('geom').feature('eastwardThinningEdges').set('entitydim', '1');
obj.m.geom('geom').feature('eastwardThinningEdges').set('add', {'eastwardThinningCumulative'});
obj.m.geom('geom').feature('eastwardThinningEdges').set('subtract', {'boundariesAdjacentToEntireSurface'});

obj.m.geom('geom').create('lateralThinningEdges', 'DifferenceSelection');
obj.m.geom('geom').feature('lateralThinningEdges').label('lateralThinningEdges');
obj.m.geom('geom').feature('lateralThinningEdges').set('entitydim', '1');
obj.m.geom('geom').feature('lateralThinningEdges').set('add', {'eastwardThinningCumulative','westwardThinningCumulative'});
obj.m.geom('geom').feature('lateralThinningEdges').set('subtract', {'eastwardThinningEdges','westwardThinningEdges'});

% meshing domain selections
% obj.m.geom('geom').selection.create('gapLeftSelection', 'CumulativeSelection');
obj.m.geom('geom').create('gapMeshingDomains', 'DifferenceSelection');
obj.m.geom('geom').feature('gapMeshingDomains').label('gapMeshingDomains');
obj.m.geom('geom').feature('gapMeshingDomains').set('entitydim', '2');
obj.m.geom('geom').feature('gapMeshingDomains').set('add', {'gapLeftSelection','gapRightSelection'});
obj.m.geom('geom').feature('gapMeshingDomains').set('subtract', {'eastwardThinningCumulative','westwardThinningCumulative'});

obj.m.geom('geom').create('domainsAdjacentToThinningAreas', 'AdjacentSelection');
obj.m.geom('geom').feature('domainsAdjacentToThinningAreas').label('domainsAdjacentToThinningAreas');
obj.m.geom('geom').feature('domainsAdjacentToThinningAreas').set('entitydim', '2');
obj.m.geom('geom').feature('domainsAdjacentToThinningAreas').set('outputdim', '2');
obj.m.geom('geom').feature('domainsAdjacentToThinningAreas').set('input', {'eastwardThinningCumulative','westwardThinningCumulative'});


obj.m.geom('geom').create('ddlDomainsAtBpeEnds', 'IntersectionSelection');
obj.m.geom('geom').feature('ddlDomainsAtBpeEnds').label('ddlDomainsAtBpeEnds');
obj.m.geom('geom').feature('ddlDomainsAtBpeEnds').set('entitydim', '2');
obj.m.geom('geom').feature('ddlDomainsAtBpeEnds').set('input', {'domainsAdjacentToThinningAreas','regionOfFirstDebyeLength'});

obj.m.geom('geom').create('ddlBoundariesAtBpeEnds', 'AdjacentSelection');
obj.m.geom('geom').feature('ddlBoundariesAtBpeEnds').label('ddlBoundariesAtBpeEnds');
obj.m.geom('geom').feature('ddlBoundariesAtBpeEnds').set('entitydim', '2');
obj.m.geom('geom').feature('ddlBoundariesAtBpeEnds').set('outputdim', '1');
obj.m.geom('geom').feature('ddlBoundariesAtBpeEnds').set('input', 'ddlDomainsAtBpeEnds');

obj.m.geom('geom').create('bpeSurfaceAtBpeEnds', 'IntersectionSelection');
obj.m.geom('geom').feature('bpeSurfaceAtBpeEnds').label('bpeSurfaceAtBpeEnds');
obj.m.geom('geom').feature('bpeSurfaceAtBpeEnds').set('entitydim', '1');
obj.m.geom('geom').feature('bpeSurfaceAtBpeEnds').set('input', {'ddlBoundariesAtBpeEnds','bpeSurfaceCumulative'});

obj.m.geom('geom').create('domainsAdjacentToEndPieces', 'AdjacentSelection');
obj.m.geom('geom').feature('domainsAdjacentToEndPieces').label('domainsAdjacentToEndPieces');
obj.m.geom('geom').feature('domainsAdjacentToEndPieces').set('entitydim', '2');
obj.m.geom('geom').feature('domainsAdjacentToEndPieces').set('outputdim', '2');
obj.m.geom('geom').feature('domainsAdjacentToEndPieces').set('input', {'gapLeftSelection','gapRightSelection'});

obj.m.geom('geom').create('upperDomainsAtBpeEnds', 'DifferenceSelection');
obj.m.geom('geom').feature('upperDomainsAtBpeEnds').label('upperDomainsAtBpeEnds');
obj.m.geom('geom').feature('upperDomainsAtBpeEnds').set('entitydim', '2');
obj.m.geom('geom').feature('upperDomainsAtBpeEnds').set('add', {'domainsAdjacentToEndPieces'});
obj.m.geom('geom').feature('upperDomainsAtBpeEnds').set('subtract', {'regionOfFirstDebyeLength'});

% domains for mapped mash at ends
obj.m.geom('geom').create('mappedMeshDomainsAtEnds', 'UnionSelection');
obj.m.geom('geom').feature('mappedMeshDomainsAtEnds').set('entitydim', '2'); % object
obj.m.geom('geom').feature('mappedMeshDomainsAtEnds').label('mappedMeshDomainsAtEnds');
obj.m.geom('geom').feature('mappedMeshDomainsAtEnds').set('input', {'ddlDomainsAtBpeEnds',...
    'westwardThinningCumulative','eastwardThinningCumulative'});
% obj.m.geom('geom').feature('domainUnion').set('contributeto', 'domainUnionCumulative');

obj.m.geom('geom').create('triangularMeshDomainsAtEnds', 'UnionSelection');
obj.m.geom('geom').feature('triangularMeshDomainsAtEnds').set('entitydim', '2'); % object
obj.m.geom('geom').feature('triangularMeshDomainsAtEnds').label('triangularMeshDomainsAtEnds');
obj.m.geom('geom').feature('triangularMeshDomainsAtEnds').set('input', {...
    'upperDomainsAtBpeEnds','gapMeshingDomains'});



obj.m.geom('geom').feature('fin').set('repairtol', '1e-10');
obj.m.geom('geom').feature('fin').set('action', 'assembly');
obj.m.geom('geom').run;
% obj.ddl = obj.m.geom('geom').feature('ddl');
% obj.m.selection('leftBoundaryOfSurface') = obj.m.geom('geom').feature('leftBoundaryOfSurface');
% obj.m.selection('rightBoundaryOfSurface') = obj.m.geom('geom').feature('rightBoundaryOfSurface');

%% selections important for extracting vertex indices
obj.leftBoundaryOfSurface = obj.m.selection.create('leftBoundaryOfSurface', 'Box');
obj.rightBoundaryOfSurface = obj.m.selection.create('rightBoundaryOfSurface', 'Box');
obj.leftBoundaryOfZetaPlane = obj.m.selection.create('leftBoundaryOfZetaPlane', 'Box');
obj.rightBoundaryOfZetaPlane = obj.m.selection.create('rightBoundaryOfZetaPlane', 'Box');

obj.m.selection('leftBoundaryOfSurface').set('entitydim', '0');
obj.m.selection('leftBoundaryOfSurface').label('leftBoundaryOfSurface');
obj.m.selection('leftBoundaryOfSurface').set('xmin', '-w_bpe/2-w_bulkLeft');
obj.m.selection('leftBoundaryOfSurface').set('xmax', '-w_bpe/2-w_bulkLeft');
obj.m.selection('leftBoundaryOfSurface').set('ymin', '0');
obj.m.selection('leftBoundaryOfSurface').set('ymax', '0');

obj.m.selection('rightBoundaryOfSurface').set('entitydim', '0');
obj.m.selection('rightBoundaryOfSurface').label('rightBoundaryOfSurface');
obj.m.selection('rightBoundaryOfSurface').set('xmin', 'w_bpe/2+w_bulkRight');
obj.m.selection('rightBoundaryOfSurface').set('xmax', 'w_bpe/2+w_bulkRight');
obj.m.selection('rightBoundaryOfSurface').set('ymin', '0');
obj.m.selection('rightBoundaryOfSurface').set('ymax', '0');

obj.m.selection('leftBoundaryOfZetaPlane').set('entitydim', '0');
obj.m.selection('leftBoundaryOfZetaPlane').label('leftBoundaryOfZetaPlane');
obj.m.selection('leftBoundaryOfZetaPlane').set('xmin', '-w_bpe/2-w_bulkLeft');
obj.m.selection('leftBoundaryOfZetaPlane').set('xmax', '-w_bpe/2-w_bulkLeft');
obj.m.selection('leftBoundaryOfZetaPlane').set('ymin', 'epsilon');
obj.m.selection('leftBoundaryOfZetaPlane').set('ymax', 'epsilon');

obj.m.selection('rightBoundaryOfZetaPlane').set('entitydim', '0');
obj.m.selection('rightBoundaryOfZetaPlane').label('rightBoundaryOfZetaPlane');
obj.m.selection('rightBoundaryOfZetaPlane').set('xmin', 'w_bpe/2+w_bulkRight');
obj.m.selection('rightBoundaryOfZetaPlane').set('xmax', 'w_bpe/2+w_bulkRight');
obj.m.selection('rightBoundaryOfZetaPlane').set('ymin', 'epsilon');
obj.m.selection('rightBoundaryOfZetaPlane').set('ymax', 'epsilon');
    
end