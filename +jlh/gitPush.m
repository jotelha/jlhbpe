function gitPush(branch)
    system('start-ssh-agent');

    if(~exist('branch','var'))
        system('git push --all origin','-echo');
    else
        system(sprintf('git push origin %s',branch),'-echo');
    end
end