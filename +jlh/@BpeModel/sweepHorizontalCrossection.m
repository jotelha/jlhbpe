function sweepHorizontalCrossection(obj,dset,stepSize,left,right,lower,upper)
    % stepSize ain fractions of L
    if ~exist('left','var')
        left = '-w_bpe/2-w_bulkLeft';
    end
    if ~exist('right','var')
        right = 'w_bpe/2+w_bulkRight';
    end
    if ~exist('lower','var')
        lower = 0;
    end
    if ~exist('upper','var')
        upper = 1;
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
        obj.plotAlongCrossection('multiPurposeCutLine',sprintf('horizontalCrossection_%.3f',pos));
    end
end