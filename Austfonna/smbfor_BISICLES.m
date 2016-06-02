
%*************************************************************************%
%**********************************Preparing SMB**************************%
%*************************************************************************%


%close all;

clear all;
delete ASF_SMB90s.nc
fprintf('Data loading...\n');
path = '/wrk/ygong/ASFforwardFiles/Invers/PostProcessing/';
filename = 'hirham_smb_1990_1999_mean.nc.txt';
smb_data=load([path,filename]);
% smb_90s = ncread([path,filename],'smb');
% 
% x_utm = ncread([path,filename],'x_utm');
% 
% y_utm = ncread([path,filename],'y_utm');

load ASF_B3_mask.mat;
%mask = Mask_ASF_Elmer;


%*************************************************************************%

% fprintf('converting 2d to colume data...\n');
% i = 0;
% for rows = 1:size(smb_90s,1)
%     smb_data(1+i*size(smb_90s,2):size(smb_90s,2)+i*size(smb_90s,2),1) = double(rot90(x_utm(rows,:),-1));
%     smb_data(1+i*size(smb_90s,2):size(smb_90s,2)+i*size(smb_90s,2),2) = double(rot90(y_utm(rows,:),-1));
%     smb_data(1+i*size(smb_90s,2):size(smb_90s,2)+i*size(smb_90s,2),3) = double(rot90(smb_90s(rows,:),-1));
%     i=i+1;
% end 

% start and end locations of 33X in UTM coordinates
xs = 545000;  xe = 755000;
ys = 8800000; ye = 8950000;
stepD = 400; % stepping through the regular array for nearest
             % neighbour calculations (in practice this gives grid
             % size in metres, you probably want somewhere between
             % 100m - 500m for Austfonna)
stepV = 1; % stepping through the scattered data (set this to 1 to
           % use all the data, or a higher integer to skip some
           % data, which will speed up this script at the cost of
           % accuracy)

fprintf('interpolating bed onto regular 2d array...\n')
[smb2d, Xq, Yq]  = interp2array(xs,xe,ys,ye,smb_data,stepV,stepD,-2000.);
figure(1),imagesc(smb2d);colorbar;

row_exp=384; colume_exp=640;
%smb_ini = flipud(smb_ini2d);
roi = 910;
smb_ini = smb2d/roi;smb_ini(row_exp,colume_exp)=0;smb_ini(mask==0.) = 0.;smb_ini(isnan(smb_ini)) = 0.;
if(isnan(smb_ini))
    fprintf('there is nans in smb_ini');
end

figure(3),imagesc(smb_ini);colorbar;hold on;
plot(ASF_index(1,:),ASF_index(2,:),'k-');
%end
if 0
fprintf('creating netcdf files...\n')
fprintf('getting the data...\n')   
y=Yq;y(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);
%x=rot90(Yq,-1);x(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);
sigma=10;

x=rot90(Xq,-1);x(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);
%y=Xq;y(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);


nccreate('ASF_SMB90s.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('ASF_SMB90s.nc','x',x);

nccreate('ASF_SMB90s.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('ASF_SMB90s.nc','y',y);

nccreate('ASF_SMB90s.nc','sigma','Dimensions',{'sigma' 10},'Format','classic');
ncwrite('ASF_SMB90s.nc','sigma',sigma);


%filename={bed_topo,'surf_topo','Vel_x','Vel_y','thickness_without_NaN','Vel_x_without_NaN','Vel_y_without_NaN'};
%nccreate('ASFdata_for_BISICLES.nc','thickness','Dimensions',{'y' 751 'x' 1051},'Format','classic');
%ncwrite('ASFdata_for_BISICLES.nc','thickness',thickness);


    data=smb_ini;filename='smb_ini';
    
    nccreate('ASF_SMB90s.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp },'Format','classic');
    ncwrite('ASF_SMB90s.nc',filename,rot90(data,-1));
    
 
ncdisp('ASF_SMB90s.nc');
end
fprintf('end\n')





