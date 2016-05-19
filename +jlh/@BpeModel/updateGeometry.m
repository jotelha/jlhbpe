function obj = updateGeometry(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    import jlh.hf.*
    import jlh.*
    
    %% modifying geometry
    fprintf('  Modifying geometry...\n');

    obj.geom.repairTol(1e-16); % should be smaller than smallest elements
    % obj.space.set('intervals','many');
    % space.set('p','0,lambdaD,L'); % an interval with three points
    % space.set('p','0,lambdaD/L,1'); % an interval with three points. dimensionless

    %electrode_surface = space.

    obj.ddl.set('size', {'w' 'epsilon'});
    obj.ddl.set('pos', {'-w/2' '0'});
%     obj.space.setIndex('layer', 'epsilon', 0);
%     obj.space.set('layername', {'dl'}); % diffuse layer
    obj.space.set('size', {'w' '1-epsilon'});
    obj.space.set('pos', {'-w/2' 'epsilon'});  
    
    obj.leftEdgeOfBpe.setIndex('p', '-w_bpe/2', 0, 0);
    obj.leftEdgeOfBpe.setIndex('p', '0', 1, 0);
    obj.rightEdgeOfBpe.setIndex('p', 'w_bpe/2', 0, 0);
    obj.rightEdgeOfBpe.setIndex('p', '0', 1, 0);
    
    if(obj.explicitElectrodeGeometry)
        obj.leftEdgeOfCathode.setIndex('p', '-w_bpe/2-w_insulatorLeft-w_cathode', 0, 0);
        obj.leftEdgeOfCathode.setIndex('p', '0', 1, 0);
        obj.rightEdgeOfCathode.setIndex('p', '-w_bpe/2-w_insulatorLeft', 0, 0);
        obj.rightEdgeOfCathode.setIndex('p', '0', 1, 0);
        
        obj.leftEdgeOfAnode.setIndex('p', 'w_bpe/2+w_insulatorRight', 0, 0);
        obj.leftEdgeOfAnode.setIndex('p', '0', 1, 0);
        obj.rightEdgeOfAnode.setIndex('p', 'w_bpe/2+w_insulatorRight+w_anode', 0, 0);
        obj.rightEdgeOfAnode.setIndex('p', '0', 1, 0);
    end
    
    obj.geom.feature('fin').set('repairtol', '1.0E-10');
    % mphgeom(model,'geom','vertexmode','on'); % plots the geometry
    %obj.geom.feature('fin').set('repairtol', '1.0E-16');
	obj.geom.run('fin');
%     obj.geom.run;
end