function sweepVerticalCrossection(obj,dset,stepSize,left,right,lower,upper)
    if ~exist('left','var')
        left = 'XleftBoundary';
    end
    if ~exist('right','var')
        right = 'XrightBoundary';
    end
    if ~exist('lower','var')
        lower = '0';
    end
    if ~exist('upper','var')
        upper = 'H';
    end
    if isnumeric(lower)
        lower = num2str(lower);
    end
    if isnumeric(upper)
        upper = num2str(upper);
    end

    obj.m.result.dataset('multiPurposeCutLine').set('data', dset);
    obj.m.result.dataset('multiPurposeCutLine').label('multiPurposeCutLine');
    obj.m.result.dataset('multiPurposeCutLine').set('spacevars', 'cx');

    obj.m.result('multiPurpose1dPlotGroup').set('data','multiPurposeCutLine');

    for pos=( left:stepSize:right )
        obj.m.result.dataset('multiPurposeCutLine').set('genpoints', { num2str(pos) lower; num2str(pos) upper});
        obj.plotAlongCrossection('multiPurposeCutLine',sprintf('verticalCrossection_%.3f',pos));
    end
end