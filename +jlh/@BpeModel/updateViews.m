function obj = updateViews(obj)
fprintf('  Update views...\n');

%     obj.standardView.set('default', 'on');
    obj.m.view('standardView').label('standardView');

    obj.m.view('standardView').set('locked', 'off');

    obj.m.view('standardView').axis.set('viewscaletype', 'manual');  
%     obj.standardView.axis.set('ymin', 0);
%     obj.standardView.axis.set('ymax', 1);
 
%     obj.standardView.set('locked', true);

    obj.m.view('standardView').axis.set('abstractviewxscale', '1');
    obj.m.view('standardView').axis.set('abstractviewyscale', obj.w/4);
    obj.m.view('standardView').axis.set('xmin', -1.5*obj.w/2);    
    obj.m.view('standardView').axis.set('xmax', 1.5*obj.w/2);
    obj.m.view('standardView').axis.set('ymin', 0);
    obj.m.view('standardView').axis.set('ymax', 1);
% %       
%     obj.standardView.axis.set('abstractviewbratio', 1);
%     obj.standardView.axis.set('abstractviewtratio', 1);
%     obj.standardView.axis.set('abstractviewrratio', 1);   
%     obj.standardView.axis.set('abstractviewlratio', 1);
    
%     obj.ddlView.set('viewscaletype','manual');

%     obj.m.view('standardView').set('yscale', obj.w/4);
    obj.m.view('standardView').set('locked', 'on');


    obj.m.view('ddlView').label('ddlView');
%     obj.ddlView.axis.set('auto', 'off'); % manual changes are not stored
%     obj.ddlView.axis.set('equal', 'off'); % use different scaling for different directions
     obj.m.view('ddlView').axis.set('viewscaletype', 'manual'); 
     obj.m.view('ddlView').axis.set('yscale', obj.w/obj.epsilon);
%      obj.ddlView.axis.set('xmin', -obj.w/2);
%      obj.ddlView.axis.set('xmax', obj.w/2);
     obj.m.view('ddlView').axis.set('ymin', 0);
     obj.m.view('ddlView').axis.set('ymax', obj.epsilon);
%      obj.ddlView.axis.set('abstractviewxscale', 1);
%      obj.ddlView.axis.set('abstractviewyscale', 1/obj.gamma);
%      obj.ddlView.axis.set('abstractviewbratio', 1);
%      obj.ddlView.axis.set('abstractviewtratio', 1/obj.gamma);
%      obj.ddlView.axis.set('abstractviewrratio', 1/2);
%      obj.ddlView.axis.set('abstractviewlratio', -1/2);
%      obj.ddlView.set('locked', true);
   
%     model.view('view3').set('locked', true);
%     model.view('view3').axis.set('abstractviewxscale', '9.332531481049955E-4');
%     model.view('view3').axis.set('ymin', '-4.893331788480282E-5');
%     model.view('view3').axis.set('yscale', '1000');
%     model.view('view3').axis.set('xmax', '-4.2489399909973145');
%     model.view('view3').axis.set('abstractviewyscale', '9.332531476502481E-7');
%     model.view('view3').axis.set('abstractviewbratio', '-4.893331788480282E-5');
%     model.view('view3').axis.set('abstractviewtratio', '-0.9995001554489136');
%     model.view('view3').axis.set('abstractviewrratio', '-0.9248939752578735');
%     model.view('view3').axis.set('xmin', '-5.2213897705078125');
%     model.view('view3').axis.set('abstractviewlratio', '-0.02213897742331028');
%     model.view('view3').axis.set('viewscaletype', 'manual');
%     model.view('view3').axis.set('ymax', '4.99819521792233E-4');
end