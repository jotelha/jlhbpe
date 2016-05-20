function gitPush(branch)
    if(~exist('branch','var'))
        system('git push --all origin');
    else
        system(sprintf('git push origin %s',branch));
    end
end