%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepera mask for VAF calculateion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete B3_mask.nc;
delete ASF_mask.nc;
clear all;
close all;

load ASF_B3_mask.mat;

B3_mask = mask;
B3_mask(mask == 1) = 0;
B3_mask(mask == 2) = 1;

ASF_mask = mask;
ASF_mask(mask == 2) = 1;

figure(1),imagesc(B3_mask),colorbar;
figure(2),imagesc(ASF_mask),colorbar;

fprintf('saving...\n');

xs = 545000;  xe = 755000;
ys = 8800000; ye = 8950000;
stepD = 400;row_exp=size(mask,1); colume_exp=size(mask,2);

Xq  = linspace(xs,xe,floor((xe-xs)/stepD+1));
Yq  = linspace(ys,ye,floor((ye-ys)/stepD+1));
Yq(end:-1:1) = Yq(1:end);

y=Yq;y(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);

sigma=1;

x=rot90(Xq,-1);x(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);
%y=Xq;y(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);

%############################
nccreate('ASF_mask.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('ASF_mask.nc','x',x);

nccreate('ASF_mask.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('ASF_mask.nc','y',y);

nccreate('ASF_mask.nc','sigma','Dimensions',{'sigma' sigma},'Format','classic');
ncwrite('ASF_mask.nc','sigma',sigma);


data=ASF_mask;filename='ASF_mask';
nccreate('ASF_mask.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp },'Format','classic');
ncwrite('ASF_mask.nc',filename,rot90(data,-1));
%############################
nccreate('B3_mask.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('B3_mask.nc','x',x);

nccreate('B3_mask.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('B3_mask.nc','y',y);

nccreate('B3_mask.nc','sigma','Dimensions',{'sigma' sigma},'Format','classic');
ncwrite('B3_mask.nc','sigma',sigma);


data=B3_mask;filename='B3_mask';
nccreate('B3_mask.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp },'Format','classic');
ncwrite('B3_mask.nc',filename,rot90(data,-1));

ncdisp('B3_mask.nc');
ncdisp('ASF_mask.nc');



