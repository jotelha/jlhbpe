function obj = setInitialValues(obj,dst_study_id,src_study_id)
%     model.study(study_id).feature('stationaryStudyStep1').set('initstudyhide', 'on');
%     model.study(study_id).feature('stationaryStudyStep1').set('initsolhide', 'on');
%     model.study(study_id).feature('stationaryStudyStep1').set('solnumhide', 'on');
%     model.study(study_id).feature('stationaryStudyStep1').set('pname', 'deltaPhiRampFactor');
%     model.study(study_id).feature('stationaryStudyStep1').set('plistarr', '1');
%     model.study(study_id).feature('stationaryStudyStep1').set('punit', '');
%     model.study(study_id).feature('stationaryStudyStep1').set('useparam', 'on');
%     model.study(study_id).feature('stationaryStudyStep1').set('pcontinuationmode', 'no');
%     model.study(study_id).feature('stationaryStudyStep1').set('sweeptype', 'filled');
    model.study(dst_study_id).feature('stationaryStudyStep1').set('geomselection', 'geom');
%     model.study(study_id).feature('stationaryStudyStep1').set('mesh', {'geom', 'standardMesh'});
    model.study(dst_study_id).feature('stationaryStudyStep1').set('useinitsol', 'on');
    model.study(dst_study_id).feature('stationaryStudyStep1').set('initmethod', 'sol');
    model.study(dst_study_id).feature('stationaryStudyStep1').set('initstudy', src_study_id);
    model.study(dst_study_id).feature('stationaryStudyStep1').set('solnum', 'auto');
end