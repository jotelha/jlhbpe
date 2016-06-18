function sweepHorizontalCrossection(obj,dset,stepSize,left,right,lower,upper,plots)
    if ~exist('left','var')
        left = 'XleftBoundary';
    end
    if ~exist('right','var')
        right = 'XrightBoundary';
    end
    if ~exist('lower','var')
        lower = 0;
    end
    if ~exist('upper','var')
        upper = obj.L;
    end
    if isnumeric(left)
        left = num2str(left);
    end
    if isnumeric(right)
        right = num2str(right);
    end
        

    obj.m.result.dataset('multiPurposeCutLine').set('data', dset);
    obj.m.result.dataset('multiPurposeCutLine').label('multiPurposeCutLine');
    obj.m.result.dataset('multiPurposeCutLine').set('spacevars', 'cx');

    obj.m.result('multiPurpose1dPlotGroup').set('data','multiPurposeCutLine');

%     for pos=( (-m.w_bpe/2-m.w_insulatorLeft):stepSize:(m.w_bpe/2+m.w_insulatorRight) )
    for pos=( lower:stepSize:upper )
        obj.m.result.dataset('multiPurposeCutLine').set('genpoints', { left num2str(pos);  right num2str(pos)});
        plotAlongCrossection('multiPurposeCutLine','x',sprintf('horizontalCrossection_%.3f',pos),plots);
    end
end