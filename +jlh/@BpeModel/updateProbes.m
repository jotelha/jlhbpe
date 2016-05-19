function updateProbes(obj,dset)

    integrationTemplate = 'integrateSurface(%s)';

    surfaceProbeTable.comments('Surface probe');

%     phiSurfaceProbe.set('data',dset);
% %     phiSurfaceProbe.selection.named('surfaceNode');
%     phiSurfaceProbe.set('probetag', 'none');
%     phiSurfaceProbe.label('potential phi surface probe');
%     phiSurfaceProbe.set('expr', 'phi');
%     phiSurfaceProbe.set('descr', 'Surface potential.');
%     phiSurfaceProbe.set('table', 'probeTable');
%     phiSurfaceProbe.set('unit', '');
%     phiSurfaceProbe.setResult;


    for i=1:numberOfSpecies  
%         cSurfaceProbe{i}.set('data',dset);
%         cSurfaceProbe{i}.selection.named('surfaceNode');
%         cSurfaceProbe{i}.set('probetag', 'none');
%         cSurfaceProbe{i}.label(sprintf('%s surface probe',c_id{i}));
%         cSurfaceProbe{i}.set('expr', c_id{i});
%         cSurfaceProbe{i}.set('descr', 'Species concentration on surface.');
%         cSurfaceProbe{i}.set('table', 'probeTable');
%         cSurfaceProbe{i}.set('unit', '');
%         cSurfaceProbe{i}.setResult;

        obj.NSurfaceProbe{i}.set('data',dset);
%         obj.NSurfaceProbe{i}.selection.named('surfaceNode');
        obj.NSurfaceProbe{i}.set('probetag', 'none');
        obj.NSurfaceProbe{i}.label(sprintf('%s surface probe',obj.N_id{i}));
        obj.NSurfaceProbe{i}.set('expr', sprintf(integrationTemplate,obj.N_id{i}));
        obj.NSurfaceProbe{i}.set('descr', 'Species flux on surface. Positive flux inwards.');
        obj.NSurfaceProbe{i}.set('table', 'probeTable');
        obj.NSurfaceProbe{i}.set('unit', '');
        obj.NSurfaceProbe{i}.setResult;
    end
    
   for i = 1:obj.nReactions
        obj.iSurfaceProbe{i}.set('data',dset);
%         obj.NSurfaceProbe{i}.selection.named('surfaceNode');
        obj.iSurfaceProbe{i}.set('probetag', 'none');
        obj.iSurfaceProbe{i}.label(sprintf('%s surface probe',obj.i_id{i}));
        obj.iSurfaceProbe{i}.set('expr', sprintf(integrationTemplate,obj.i_id{i}));
        obj.iSurfaceProbe{i}.set('descr', 'Current due to single reaction');
        obj.iSurfaceProbe{i}.set('table', 'probeTable');
        obj.iSurfaceProbe{i}.set('unit', '');
        obj.iSurfaceProbe{i}.setResult;
    end

    iCathodicSurfaceProbe.set('data',dset);
%     iCathodicSurfaceProbe.selection.named('surfaceNode');
    iCathodicSurfaceProbe.set('probetag', 'none');
    iCathodicSurfaceProbe.label('i_cathodic surface probe');
    iCathodicSurfaceProbe.set('expr', sprintf(integrationTemplate,'i_cathodic'));
    iCathodicSurfaceProbe.set('descr', 'Cathodic current due to surface reactions. A positive current means negative charge flux into the electrolyte domain.');
    iCathodicSurfaceProbe.set('table', 'probeTable');
    iCathodicSurfaceProbe.set('unit', '');
    iCathodicSurfaceProbe.setResult;

    iAnodicSurfaceProbe.set('data',dset);
%     iAnodicSurfaceProbe.selection.named('surfaceNode');
    iAnodicSurfaceProbe.set('probetag', 'none');
    iAnodicSurfaceProbe.label('i_anodic surface probe');
    iAnodicSurfaceProbe.set('expr', sprintf(integrationTemplate,'i_anodic'));
    iAnodicSurfaceProbe.set('descr', 'Anodic current due to surface reactions. A positive current means negative charge flux into the electrolyte domain.');
    iAnodicSurfaceProbe.set('table', 'probeTable');
    iAnodicSurfaceProbe.set('unit', '');
    iAnodicSurfaceProbe.setResult;

    iTotalSurfaceProbe.set('data',dset);
%     iTotalSurfaceProbe.selection.named('surfaceNode');
    iTotalSurfaceProbe.set('probetag', 'none');
    iTotalSurfaceProbe.label('i_total surface probe');
    iTotalSurfaceProbe.set('expr', sprintf(integrationTemplate,'i_total'));
    iTotalSurfaceProbe.set('descr', 'Total current due to surface reactions. A positive current means negative charge flux into the electrolyte domain.');
    iTotalSurfaceProbe.set('table', 'probeTable');
    iTotalSurfaceProbe.set('unit', '');
    iTotalSurfaceProbe.setResult;
end
end