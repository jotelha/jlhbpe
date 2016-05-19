function obj = initializeExporter(obj)
    obj.m.result.export('plotExporter1d').set('printunit', 'mm');
    obj.m.result.export('plotExporter1d').set('webunit', 'px');
    obj.m.result.export('plotExporter1d').set('printheight', '90');
    obj.m.result.export('plotExporter1d').set('webheight', '600');
    obj.m.result.export('plotExporter1d').set('printwidth', '120');
    obj.m.result.export('plotExporter1d').set('webwidth', '800');
    obj.m.result.export('plotExporter1d').set('printlockratio', 'off'); %  off
    obj.m.result.export('plotExporter1d').set('weblockratio', 'off'); % off
    obj.m.result.export('plotExporter1d').set('printresolution', '300');
    obj.m.result.export('plotExporter1d').set('webresolution', '96');
    obj.m.result.export('plotExporter1d').set('size', 'manualprint'); % current, manualprint, manualweb
    obj.m.result.export('plotExporter1d').set('antialias', 'on');
    obj.m.result.export('plotExporter1d').set('zoomextents', 'on'); % off
    obj.m.result.export('plotExporter1d').set('title', 'off');
%     obj.m.result.export('plotExporter1d') .set('legend', 'on');
    obj.m.result.export('plotExporter1d').set('legend', 'off');
    obj.m.result.export('plotExporter1d').set('logo', 'off');
    obj.m.result.export('plotExporter1d').set('options', 'on');
    obj.m.result.export('plotExporter1d').set('fontsize', '9');
    obj.m.result.export('plotExporter1d').set('customcolor', [1 1 1]);
    obj.m.result.export('plotExporter1d').set('background', 'color');
    obj.m.result.export('plotExporter1d').set('axes', 'on');
%     obj.m.result.export('plotExporter1d') .set('grid', 'on'); % jlh, for 3d
    obj.m.result.export('plotExporter1d').set('qualitylevel', '92');
    obj.m.result.export('plotExporter1d').set('qualityactive', 'on'); % off
    obj.m.result.export('plotExporter1d').set('imagetype', 'png');
end