function nccreat4BISICLES(xq,yq,sigma,ncfilename,structuredata)
if ~(isstruct((structuredata)))
    fprintf('the data all needs to be saved in one structure data with names,such as structuredata.data, structuredata.name! \n');
    return;
end
%row_exp=384; colume_exp=640;

y=yq;%y(end:row_exp) = yq(end,end):400:(yq(end,end)+(row_exp-size(yq,2))*400);



x=rot90(xq,-1);%x(end:colume_exp) = xq(end,end):400:(xq(end,end)+(colume_exp-size(xq,2))*400);



nccreate(ncfilename,'x','Dimensions',{'x' length(x) },'Format','classic');
ncwrite(ncfilename,'x',x);

nccreate(ncfilename,'y','Dimensions',{'y' length(y)},'Format','classic');
ncwrite(ncfilename,'y',y);

nccreate(ncfilename,'sigma','Dimensions',{'sigma' 10},'Format','classic');
ncwrite(ncfilename,'sigma',sigma);
for i = 1:length(structuredata)

data = structuredata(i).data;
filename=structuredata(i).name;
nccreate(ncfilename,filename,'Dimensions',{'x' length(x) 'y' length(y) },'Format','classic');
ncwrite(ncfilename,filename,rot90(data,-1));
end
fprintf(['file ', ncfilename, ' created.\n']);
end
