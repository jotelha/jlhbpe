function obj = runSampleProject(obj)
    import com.comsol.model.*
    import com.comsol.model.util.*
    obj.zeroSurfaceCurrent.active(false);
    ModelUtil.ShowProgress = true;
    obj.stationaryStudy.run();
    ModelUtil.ShowProgress = false;
end