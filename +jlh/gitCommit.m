function gitCommit(msg)
    system('git add *.m','-echo');
    %system('git add img/*','-echo');
    %system('git add dat/*','-echo');
    %system('git add TODO','-echo');
    system(['git commit -a -m "',msg,'"'],'-echo');
end