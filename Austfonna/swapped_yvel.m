close all;
clear all;
load('ASFdata_2011Ori_for_BISICLES.mat');
load vvel_asf.dat;
load uvel_asf.dat;

x=Xq;sigma=10;
for row= 1:length(Yq)
    y(row,1)=Yq(1,row);
end

data=-vvel_asf;filename='yvel_neg';
row_exp=384; colume_exp=640;

figure(1),imagesc(vvel_asf),colorbar;
figure(2),imagesc(data),quiver(uvel_asf,vvel_asf),colorbar;
%figure(3),quiver(uvel_asf,vvel_asf

nccreate('ASF_yvel_neg.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('ASF_yvel_neg.nc','x',x);

nccreate('ASF_yvel_neg.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('ASF_yvel_neg.nc','y',y);

nccreate('ASF_yvel_neg.nc','sigma','Dimensions',{'sigma' 10},'Format','classic');
ncwrite('ASF_yvel_neg.nc','sigma',sigma);


nccreate('ASF_yvel_neg.nc',filename,'Dimensions',{'y' row_exp 'x' colume_exp},'Format','classic');
ncwrite('ASF_yvel_neg.nc',filename,data);

ncdisp('ASF_yvel_neg.nc');

fprintf('end\n')