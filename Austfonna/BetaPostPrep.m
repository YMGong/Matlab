%*************************************************************************%
% smb data is from the merged data created by Rupert's m scrtipt, not the %
% original data!!                                                         %
%The data must have utm coordinate infor!!                                %
%*************************************************************************%
close all;
clear all;

fprintf('loading data...\n')
pathB='/wrk/ygong/BISICLES/ASF_data/';
fnameB='ASF_2011_vel_Control.mat';
pathD= '/wrk/ygong/BISICLES/ASF_data/';
fnameD='ASFdata_2011Ori_for_BISICLES.mat';
load([pathB,fnameB]);load([pathD,fnameD]);

data=Bate_2011;filename='Bate_2011';
imagesc(data);colorbar;

row_exp=384; colume_exp=640; 

fprintf('creating netcdf files...\n')
fprintf('getting the data...\n')   
y=Yq;y(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);
%x=rot90(Yq,-1);x(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);
sigma=10;

x=rot90(Xq,-1);x(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);
%y=Xq;y(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);


nccreate('ASFBate_2011.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('ASFBate_2011.nc','x',x);

nccreate('ASFBate_2011.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('ASFBate_2011.nc','y',y);

nccreate('ASFBate_2011.nc','sigma','Dimensions',{'sigma' 10},'Format','classic');
ncwrite('ASFBate_2011.nc','sigma',sigma);

nccreate('ASFBate_2011.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp},'Format','classic');
ncwrite('ASFBate_2011.nc',filename,data);
    
ncdisp('ASFBate_2011.nc');

fprintf('end\n')











