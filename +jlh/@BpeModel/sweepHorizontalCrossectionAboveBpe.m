function sweepHorizontalCrossectionAboveBpe(obj,dset,stepSize)
    % stepSize ain fractions of L

    obj.m.result.dataset('multiPurposeCutLine').set('data', dset);
    obj.m.result.dataset('multiPurposeCutLine').label('multiPurposeCutLine');
    obj.m.result.dataset('multiPurposeCutLine').set('spacevars', 'cx');

    obj.m.result('multiPurpose1dPlotGroup').set('data','multiPurposeCutLine');

%     for pos=( (-m.w_bpe/2-m.w_insulatorLeft):stepSize:(m.w_bpe/2+m.w_insulatorRight) )
    for pos=( 0:stepSize:1 )
        obj.m.result.dataset('multiPurposeCutLine').set('genpoints', { '-w_bpe/2' num2str(pos); 'w_bpe/2' num2str(pos)});
        obj.plotAlongCrossection('multiPurposeCutLine',sprintf('horizontalCrossectionAboveBpe_%.3f',pos));
    end
end