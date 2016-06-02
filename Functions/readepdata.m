
%-------------- Reading elmerpost files by line---------------------%
% l_s is the starting line; l_e is the ending line; 
%colume_n is total number of the colume of the data in the .ep file
%specify colume_s when the starting colume is a string
%-------------------------------------------------------------------%
function data = readepdata(filename,l_s,l_e,colume_n,colume_s)

 fid=fopen(filename, 'r');
 if fid == -1;
    fprintf('\nERROR : The file do not exist!\n');
    return;
 end
 
 % ignore the lines before the starting line
 
 for i = 1:(l_s-1)
     ln(i).ignore = fgetl(fid);
 end
 
 % get the data wanted
 if colume_s == 0
 for j = l_s:l_e
     data((j-l_s+1),:)=fscanf(fid,'%f',colume_n);
 end
 else 
     for j = l_s:l_e
     name((j-l_s+1),:)=fscanf(fid,'%s',colume_s);
     data((j-l_s+1),:)=fscanf(fid,'%f',colume_n);
     end
 end
fclose(fid);
end