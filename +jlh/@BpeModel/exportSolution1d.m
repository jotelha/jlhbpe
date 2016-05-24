function txtFileName = exportSolution1d(obj,dset,dsetFolder)
    if ~exist('dsetFolder','var')
        dsetFolder = dset;
    end

    import jlh.*
    import jlh.hf.*

    fullProjectPath = [pwd(),'\',obj.projectPath];


    expressions = [{'phi'}, obj.c_id];

        
    solution1dSubFolder = [fullProjectPath,'\solution1d'];
    if( ~exist( solution1dSubFolder,'dir') )
        mkdir(solution1dSubFolder);
    end

    dsetSubFolder = [solution1dSubFolder,'\',dsetFolder];
    if( ~exist( dsetSubFolder,'dir') )
        mkdir(dsetSubFolder);
    end


%     csvFileNameTemplate = [strrep(dsetSubFolder,'\','\\'),'\\%s.csv'];
            
%     csvFileName = sprintf(csvFileNameTemplate,selections{m});
     txtFileName = [dsetSubFolder,'\solution1d.txt'];

     e = ['x',flattenCell(expressions)];
 

    d = mpheval(obj.m,e,'dataset',dset,'dataonly','on');
    v = arrayfun(@(f) (d.(sprintf('d%d',f)))', 1:numel(e), 'UniformOutput', false);

    out = cell2mat(v);
%         fid=fopen(csvFileName,'w');           
%         fprintf(fid,'%s,',e{1:(end-1)});
%         fprintf(fid,'%s\n',e{end});
%         fclose(fid);
%         dlmwrite(txtFileName,out,'-append','delimiter',',','precision',8);
    dlmwrite(txtFileName,out,'delimiter','\t','precision','%.8e');
end
              
                    