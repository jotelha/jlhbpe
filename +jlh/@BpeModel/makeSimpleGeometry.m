function obj = makeSimpleGeometry(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% modifying geometry

    obj.geom = obj.m.geom.create('geom',2); % 1 for one dimension

    obj.m.geom('geom').repairTol(1e-10); % should be smaller than smallest elements

    obj.ddl = obj.geom.create('ddl', 'Rectangle');
    obj.m.geom('geom').feature('ddl').label('ddl');
    obj.m.geom('geom').feature('ddl').set('size', {'w' 'epsilon'});
    obj.m.geom('geom').feature('ddl').set('pos', {'-w_bpe/2-w_bulkLeft' '0'});
    obj.m.geom('geom').feature('ddl').set('selresult', 'on');
    obj.m.geom('geom').feature('ddl').set('selresultshow', 'all');


    %     obj.space.setIndex('layer', 'epsilon', 0);
    %     obj.space.set('layername', {'dl'}); % diffuse layer
    obj.space = obj.geom.create('space', 'Rectangle');  
    obj.m.geom('geom').feature('space').label('space');    
    obj.m.geom('geom').feature('space').set('size', {'w' '1-epsilon'});
    obj.m.geom('geom').feature('space').set('pos', {'-w_bpe/2-w_bulkLeft' 'epsilon'});  
    obj.m.geom('geom').feature('space').set('selresult', 'on');
    obj.m.geom('geom').feature('space').set('selresultshow', 'all');

    %     obj.leftEdgeOfBpe.setIndex('p', '-w_bpe/2', 0, 0);
    %     obj.leftEdgeOfBpe.setIndex('p', '0', 1, 0);
    %     obj.rightEdgeOfBpe.setIndex('p', 'w_bpe/2', 0, 0);
    %     obj.rightEdgeOfBpe.setIndex('p', '0', 1, 0);

%     obj.m.geom('geom').feature('fin').set('repairtol', '1.0E-10');
        % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
        %obj.geom.feature('fin').set('repairtol', '1.0E-16');
        obj.geom.run('fin');
    %     obj.geom.run;

    obj.m.geom('geom').repairTol(1.0E-10);

%% cumulative selections
%     obj.m.geom('geom').selection.create('upperBoundaryCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('upperBoundaryCumulative').label('upperBoundaryCumulative');
%     obj.m.geom('geom').selection.create('leftBoundarySelection', 'CumulativeSelection');
%     obj.m.geom('geom').selection('leftBoundarySelection').label('leftBoundarySelection');
%     obj.m.geom('geom').selection.create('meshChopPrototype', 'CumulativeSelection');
%     obj.m.geom('geom').selection('meshChopPrototype').label('meshChopPrototype');
%     obj.m.geom('geom').selection.create('rightBoundarySelection', 'CumulativeSelection');
%     obj.m.geom('geom').selection('rightBoundarySelection').label('rightBoundarySelection');
%     obj.m.geom('geom').selection.create('gapLeftSelection', 'CumulativeSelection');
%     obj.m.geom('geom').selection('gapLeftSelection').label('gapLeftSelection');
%     obj.m.geom('geom').selection.create('gapRightSelection', 'CumulativeSelection');
%     obj.m.geom('geom').selection('gapRightSelection').label('gapRightSelection');
%     obj.m.geom('geom').selection.create('zetaPlaneCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('zetaPlaneCumulative').label('zetaPlaneCumulative');
%     obj.m.geom('geom').selection.create('bpeSurfaceCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('bpeSurfaceCumulative').label('bpeSurfaceCumulative');
%     obj.m.geom('geom').selection.create('entireSurfaceCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('entireSurfaceCumulative').label('entireSurfaceCumulative');
%     obj.m.geom('geom').selection.create('regionOfFirstDebyeLength', 'CumulativeSelection');
%     obj.m.geom('geom').selection('regionOfFirstDebyeLength').label('regionOfFirstDebyeLength');
% 
%     obj.m.geom('geom').selection.create('meshChopArrayCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('meshChopArrayCumulative').label('meshChopArrayCumulative');
%     obj.m.geom('geom').selection.create('domainUnionCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('domainUnionCumulative').label('domainUnionCumulative');
% 
%     obj.m.geom('geom').selection.create('unitCellPrototypeCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('unitCellPrototypeCumulative').label('unitCellPrototypeCumulative');
%     obj.m.geom('geom').selection.create('unitCellPrototypeBoundaryCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('unitCellPrototypeBoundaryCumulative').label('unitCellPrototypeBoundaryCumulative');
% 
%     obj.m.geom('geom').selection.create('eastwardThinningCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('eastwardThinningCumulative').label('eastwardThinningCumulative');
%     obj.m.geom('geom').selection.create('westwardThinningCumulative', 'CumulativeSelection');
%     obj.m.geom('geom').selection('westwardThinningCumulative').label('westwardThinningCumulative');
% 


% points
% obj.m.geom('geom').create('leftBoundaryOfZetaPlanePoint', 'Point');
% obj.m.geom('geom').feature('leftBoundaryOfZetaPlanePoint').label('leftBoundaryOfZetaPlanePoint');
% obj.m.geom('geom').feature('leftBoundaryOfZetaPlanePoint').setIndex('p', '-w_bpe/2-w_bulkLeft', 0, 0);
% obj.m.geom('geom').feature('leftBoundaryOfZetaPlanePoint').setIndex('p', 'epsilon', 1, 0);
% obj.leftBoundaryOfZetaPlane = obj.m.geom('geom').create('leftBoundaryOfZetaPlane', 'ExplicitSelection');
% obj.m.geom('geom').feature('leftBoundaryOfZetaPlane').label('leftBoundaryOfZetaPlane');
% obj.m.geom('geom').feature('leftBoundaryOfZetaPlane').selection('selection').init(0);
% obj.m.geom('geom').feature('leftBoundaryOfZetaPlane').selection('selection').set('leftBoundaryOfZetaPlanePoint(1)', [1]);
% obj.m.geom('geom').create('rightBoundaryOfZetaPlanePoint', 'Point');
% obj.m.geom('geom').feature('rightBoundaryOfZetaPlanePoint').label('rightBoundaryOfZetaPlanePoint');
% obj.m.geom('geom').feature('rightBoundaryOfZetaPlanePoint').setIndex('p', 'w_bpe/2+w_bulkRight', 0, 0);
% obj.m.geom('geom').feature('rightBoundaryOfZetaPlanePoint').setIndex('p', 'epsilon', 1, 0);
% obj.rightBoundaryOfZetaPlane = obj.m.geom('geom').create('rightBoundaryOfZetaPlane', 'ExplicitSelection');
% obj.m.geom('geom').feature('rightBoundaryOfZetaPlane').label('rightBoundaryOfZetaPlane');
% obj.m.geom('geom').feature('rightBoundaryOfZetaPlane').selection('selection').init(0);
% obj.m.geom('geom').feature('rightBoundaryOfZetaPlane').selection('selection').set('rightBoundaryOfZetaPlanePoint(1)', [1]);

% selections
obj.m.geom('geom').create('entireSurface', 'BoxSelection');
obj.m.geom('geom').feature('entireSurface').label('entireSurface');
obj.m.geom('geom').feature('entireSurface').set('entitydim', '1');
% obj.m.geom('geom').feature('entireSurface').set('ymin', '-epsilon/2');
obj.m.geom('geom').feature('entireSurface').set('ymax', 'epsilon/2');
obj.m.geom('geom').feature('entireSurface').set('condition', 'inside');

obj.m.geom('geom').create('upperBoundary', 'BoxSelection');
obj.m.geom('geom').feature('upperBoundary').label('upperBoundarySelection');
obj.m.geom('geom').feature('upperBoundary').set('entitydim', '1');
obj.m.geom('geom').feature('upperBoundary').set('ymin', '1-epsilon/2');
obj.m.geom('geom').feature('upperBoundary').set('condition', 'inside');

obj.m.geom('geom').create('leftBoundarySelection', 'BoxSelection');
obj.m.geom('geom').feature('leftBoundarySelection').label('leftBoundarySelection');
obj.m.geom('geom').feature('leftBoundarySelection').set('entitydim', '1');
obj.m.geom('geom').feature('leftBoundarySelection').set('xmax', '-w_bpe/2-w_bulkLeft+epsilon/2');
obj.m.geom('geom').feature('leftBoundarySelection').set('condition', 'inside');

obj.m.geom('geom').create('rightBoundarySelection', 'BoxSelection');
obj.m.geom('geom').feature('rightBoundarySelection').label('rightBoundarySelection');
obj.m.geom('geom').feature('rightBoundarySelection').set('entitydim', '1');
obj.m.geom('geom').feature('rightBoundarySelection').set('xmin', 'w_bpe/2+w_bulkLeft-epsilon/2');
obj.m.geom('geom').feature('rightBoundarySelection').set('condition', 'inside');

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
obj.m.geom('geom').feature('bulkBoundary').set('input', {'rightBoundarySelection' 'leftBoundarySelection' 'upperBoundary'});

obj.bpeSurface = obj.m.geom('geom').create('bpeSurface', 'UnionSelection');
obj.m.geom('geom').feature('bpeSurface').set('entitydim', '1');
obj.m.geom('geom').feature('bpeSurface').label('bpeSurface');
obj.m.geom('geom').feature('bpeSurface').set('input', {'entireSurface'});

obj.m.geom('geom').create('zetaPlaneDegenerate', 'BoxSelection');
obj.m.geom('geom').feature('zetaPlaneDegenerate').label('zetaPlaneDegenerate');
obj.m.geom('geom').feature('zetaPlaneDegenerate').set('entitydim', '1');
obj.m.geom('geom').feature('zetaPlaneDegenerate').set('ymin', 'epsilon/2');
obj.m.geom('geom').feature('zetaPlaneDegenerate').set('ymax', '3*epsilon/2');
obj.m.geom('geom').feature('zetaPlaneDegenerate').set('condition', 'inside');

obj.m.geom('geom').create('zetaPlane', 'DifferenceSelection');
obj.m.geom('geom').feature('zetaPlane').label('zetaPlane');
obj.m.geom('geom').feature('zetaPlane').set('entitydim', '1');
obj.m.geom('geom').feature('zetaPlane').set('add', 'zetaPlaneDegenerate');
obj.m.geom('geom').feature('zetaPlane').set('subtract', 'space');

% obj.zetaPlane = obj.m.geom('geom').create('zetaPlane', 'UnionSelection');
% obj.m.geom('geom').feature('zetaPlane').set('entitydim', '1');
% obj.m.geom('geom').feature('zetaPlane').label('zetaPlane');
% obj.m.geom('geom').feature('zetaPlane').set('input', {'zetaPlaneCumulative'});
% obj.entireSurface = obj.m.geom('geom').create('entireSurface', 'UnionSelection');
% obj.m.geom('geom').feature('entireSurface').set('entitydim', '1');
% obj.m.geom('geom').feature('entireSurface').label('entireSurface');
% obj.m.geom('geom').feature('entireSurface').set('input', {'entireSurfaceCumulative'});
% obj.upperBoundary = obj.m.geom('geom').create('upperBoundary', 'UnionSelection');
% obj.m.geom('geom').feature('upperBoundary').set('entitydim', '1');
% obj.m.geom('geom').feature('upperBoundary').label('upperBoundary');
% obj.m.geom('geom').feature('upperBoundary').set('input', {'upperBoundaryCumulative'});

%% selections important for unit cell meshing

% obj.m.geom('geom').create('lateralThinningEdges', 'DifferenceSelection');
% obj.m.geom('geom').feature('lateralThinningEdges').label('lateralThinningEdges');
% obj.m.geom('geom').feature('lateralThinningEdges').set('entitydim', '1');
% obj.m.geom('geom').feature('lateralThinningEdges').set('add', {'lateralBoundary'});
% obj.m.geom('geom').feature('lateralThinningEdges').set('subtract', {'eastwardThinningEdges','westwardThinningEdges', 'allExceptDdl'});


obj.m.geom('geom').create('lateralThinningEdges', 'BoxSelection');
obj.m.geom('geom').feature('lateralThinningEdges').label('lateralThinningEdges');
obj.m.geom('geom').feature('lateralThinningEdges').set('entitydim', '1');
obj.m.geom('geom').feature('lateralThinningEdges').set('ymin', 'epsilon/2');
obj.m.geom('geom').feature('lateralThinningEdges').set('ymax', 'epsilon/2');
obj.m.geom('geom').feature('lateralThinningEdges').set('condition', 'intersects');

obj.m.geom('geom').feature('fin').set('repairtol', '1e-10');
obj.m.geom('geom').feature('fin').set('action', 'assembly');
obj.m.geom('geom').run;
% obj.ddl = obj.m.geom('geom').feature('ddl');
% obj.m.selection('leftBoundaryOfSurface') = obj.m.geom('geom').feature('leftBoundaryOfSurface');
% obj.m.selection('rightBoundaryOfSurface') = obj.m.geom('geom').feature('rightBoundaryOfSurface');

%% selections important for extracting vertex indices
obj.leftBoundaryOfSurface = obj.m.selection.create('leftBoundaryOfSurface', 'Box');
obj.rightBoundaryOfSurface = obj.m.selection.create('rightBoundaryOfSurface', 'Box');
obj.leftBoundaryOfZetaPlane = obj.m.selection.create('leftBoundaryOfZetaPlane', 'Difference');
obj.rightBoundaryOfZetaPlane = obj.m.selection.create('rightBoundaryOfZetaPlane', 'Difference');
obj.m.selection.create('leftBoundaryOfZetaPlaneDegenerate', 'Box');
obj.m.selection.create('rightBoundaryOfZetaPlaneDegenerate', 'Box');

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

obj.m.selection('leftBoundaryOfZetaPlaneDegenerate').set('entitydim', '0');
obj.m.selection('leftBoundaryOfZetaPlaneDegenerate').label('leftBoundaryOfZetaPlaneDegenerate');
obj.m.selection('leftBoundaryOfZetaPlaneDegenerate').set('xmin', '-w_bpe/2-w_bulkLeft');
obj.m.selection('leftBoundaryOfZetaPlaneDegenerate').set('xmax', '-w_bpe/2-w_bulkLeft');
obj.m.selection('leftBoundaryOfZetaPlaneDegenerate').set('ymin', 'epsilon');
obj.m.selection('leftBoundaryOfZetaPlaneDegenerate').set('ymax', 'epsilon');

obj.m.selection('leftBoundaryOfZetaPlane').label('leftBoundaryOfZetaPlane');
obj.m.selection('leftBoundaryOfZetaPlane').set('entitydim', '0');
obj.m.selection('leftBoundaryOfZetaPlane').set('add', 'leftBoundaryOfZetaPlaneDegenerate');
obj.m.selection('leftBoundaryOfZetaPlane').set('subtract', 'geom_space_pnt');

obj.m.selection('rightBoundaryOfZetaPlaneDegenerate').set('entitydim', '0');
obj.m.selection('rightBoundaryOfZetaPlaneDegenerate').label('rightBoundaryOfZetaPlaneDegenerate');
obj.m.selection('rightBoundaryOfZetaPlaneDegenerate').set('xmin', 'w_bpe/2+w_bulkRight');
obj.m.selection('rightBoundaryOfZetaPlaneDegenerate').set('xmax', 'w_bpe/2+w_bulkRight');
obj.m.selection('rightBoundaryOfZetaPlaneDegenerate').set('ymin', 'epsilon');
obj.m.selection('rightBoundaryOfZetaPlaneDegenerate').set('ymax', 'epsilon');

obj.m.selection('rightBoundaryOfZetaPlane').label('rightBoundaryOfZetaPlane');
obj.m.selection('rightBoundaryOfZetaPlane').set('entitydim', '0');
obj.m.selection('rightBoundaryOfZetaPlane').set('add', 'rightBoundaryOfZetaPlaneDegenerate');
obj.m.selection('rightBoundaryOfZetaPlane').set('subtract', 'geom_space_pnt');


%% swap identity pair. source: fine mesh, destination: coarse mesh
% ap = obj.getIdentityPairsForComponent('comp1');
% for i = 1:numel(ap)
%     obj.m.pair(ap{i}).swap();
% end


    
end