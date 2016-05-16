
%*************************************************************************%
%****************** Comparision between Elmer/Ice and BISICLES ************
%*************************************************************************%

close all;
clear all;

fprintf('Loading data... \n');
%addpath /wrk/ygong/ASFforwardFiles/Invers/postprocessing/E_results/
%addpath /wrk/ygong/BISICLES/ASF_simulation/PostProcessing/B_resutlts
path2 = 'J:\ASF_data\B_results\';
%addpath /wrk/ygong/BISICLES/ASF_simulation/PostProcessing/
% load('ASFForward_ini_Elmer.mat');
% load('ASFForward_Elmer');
% load('ASFForward_BISICLES.mat');
% load('ASFForward_ini_BISICLES.mat');
load ASF_B3_mask.mat;
load profile_B3.mat;
proN = load('profile_B3_N.mat');
proM = load('profile_B3_M.mat');
proS = load('profile_B3_S.mat');
%*************************************************************************%

fprintf('Preparing data... \n');
obvel2011x = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'xvel');
obvel2011x = rot90(obvel2011x);
obvel2011y = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'yvel');
obvel2011y = rot90(obvel2011y);
velo2011 = sqrt(obvel2011x.^2 + obvel2011y.^2);
velo2011(mask == 0)= NaN;
velo2011(velo2011 == 0) = NaN;

base_name_Elmer = 'Elmer';

%volume_name = {'asf','b3'};
time_name = {'95to11'};
%smb_name = {'elevcorr'};
%E_smb_name = {'elecorr'};
%beta_name = {'celmer_stepc95','celmer_stepc2011','cbisicles_linear'};
%time_name = {'12to13'};
smb_name = {'elevcorr.ep10'};%,'noelevcorr','20thsmb'};
%smb_name = {'20thsmb_yr1301.ep6'};
%smb_name = {'20thsmb.ep6'};
%E_smb_name = {'elecorr'};
%beta_name = {'celmer_triangelnofrc','celmer_triangellargenofrc'};%,'celmer_stepc95','celmer_stepc2011'};
%beta_name = {'celmer_triangelnofrc','celmer_triangellargenofrc'};
%beta_name = {'celmer_extropln'};
beta_name = {'celmer_linear'};
temp_name = {'temp95'};
B_beta_ini_name = {'20yr'};
E_beta_ini_name = {'nothermal_ebet'};
B_path_name = 'J:\ASF_data\B_results';
E_path_name = 'J:\ASF_data\E_results';
load 'J:\ASF_data\asf_b3_mask.txt';
%*****************************configure***********************************%
%*****************************for triangles***********************************%
row_exp=384; colume_exp=640;
%figure,imagesc(mask),colorbar;
mask_tmp = mask;
mask_tmp(mask==1)=-1999;
mask_tmp(mask==2)=-2999;
%mask(252:316,397:450)=-2999;
%figure,imagesc(mask),colorbar;
B3_mask_tmp = asf_b3_mask(:,1:2);
%B3_mask_tmp(:,3) = reshape(mask_tmp,size(mask_tmp,1)*size(mask_tmp,2),1);
%dlmwrite('asf_b3_mask_extend.txt',B3_mask_tmp,'delimiter',' ','newline','pc','precision','%8.6f');
B3_mask_coor_x = reshape(asf_b3_mask(:,1),row_exp,colume_exp);
B3_mask_coor_y = reshape(asf_b3_mask(:,2),row_exp,colume_exp);
%find the index
I_xllarge = find(B3_mask_coor_x(1,:)==703400);%mask large & 2
I_xrlarge = find(B3_mask_coor_x(1,:)==723000);%mask large & 2
I_yularge = find(B3_mask_coor_y(:,1)==8840000); %mask large
I_ydlarge = find(B3_mask_coor_y(:,1)==8825200);%mask large & 2

I_xl = find(B3_mask_coor_x(1,:)==705400);%mask 
I_xr = find(B3_mask_coor_x(1,:)==723000);%mask
%I_yu = find(B3_mask_coor_y(:,1)==8840000); %mask 
I_yu = find(B3_mask_coor_y(:,1)==8834000); %mask
I_yd = find(B3_mask_coor_y(:,1)==8825200);%mask


%xaxi = [150 500];yaxi = [50 370];
xaxi = [305 445];yaxi = [180 320];
xlab = (15:20:135)*0.400; ylab = (140:-20:0)*0.400;             
gridon = 1; lineon = 0;boxon = 1;
axi_units = 'km';
%colorbarlim=[0,5000;0,15000;-10,750];
colorbarlim = [0,100;0,40];
colorbar_lab = (0:1)*100;
grey = [0.6 0.6 0.6];
%c_units = 'velocity (m day^{-1})';
%c_units = 'm a^{-1}';
c_units = '%';
%*****************************BISICLES************************************%
count = 1;
for i = 1: length(beta_name)
for j = 1: length(smb_name)

%*******************************Elmer**************************************%
filename = [E_path_name,'\',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,cell2mat(E_beta_ini_name(i)),'.ep.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_','20thsmb.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'_','ini.ep9.mat'];
%filename_ini = [E_path_name,'/','Elmer_95to11_temp95_celmer_linear_elevcorr.ep9.mat'];

load(cell2mat(filename));
%Eini=load(cell2mat(filename_ini));
%Eini=load(filename_ini);smb_name = {'elevcorr','noelevcorr'};smb_name = {'elevcorr','noelevcorr'};
Elmer_xvel_bed = matrixdata(2).data_layer(1).data;
Elmer_yvel_bed = matrixdata(3).data_layer(1).data;
Elmer_zvel_bed = matrixdata(4).data_layer(1).data;
Elmer_xvel_surf = matrixdata(2).data_layer(2).data;
Elmer_yvel_surf = matrixdata(3).data_layer(2).data;
Elmer_zvel_surf = matrixdata(4).data_layer(2).data;
Elmer_height_surf = matrixdata(1).data_layer(1).data;
Elmer_beta_surf = matrixdata(6).data_layer(1).data;

% Elmer_xvel_ini_bed = Eini.matrixdata(2).data_layer(1).data;
% Elmer_yvel_ini_bed = Eini.matrixdata(3).data_layer(1).data;
% Elmer_zvel_ini_bed = Eini.matrixdata(4).data_layer(1).data;
% Elmer_xvel_ini_surf = Eini.matrixdata(2).data_layer(2).data;
% Elmer_yvel_ini_surf = Eini.matrixdata(3).data_layer(2).data;
% Elmer_zvel_ini_surf = Eini.matrixdata(4).data_layer(2).data;
% Elmer_height_ini_surf = Eini.matrixdata(1).data_layer(1).data;
% Elmer_beta_ini_surf = Eini.matrixdata(8).data_layer(1).data;

%Elmer_height_ini_surf (mask==0) = NaN;
Elmer_height_surf (mask==0) = NaN;
%Elmer_thickness_change = Elmer_height_surf - Elmer_height_ini_surf;

bVel_Elmer = sqrt(Elmer_xvel_bed.^2 +Elmer_yvel_bed.^2+Elmer_zvel_bed.^2);
bVel_Elmer (mask==0) = NaN;

%bVel_ini_Elmer = sqrt(Elmer_xvel_ini_bed.^2 + Elmer_yvel_ini_bed.^2+Elmer_zvel_ini_bed.^2);
%bVel_ini_Elmer (mask==0) = NaN;
%Elmer_bVel_change = bVel_Elmer - bVel_ini_Elmer;

fVel_Elmer = sqrt(Elmer_xvel_surf.^2 +Elmer_yvel_surf.^2+Elmer_zvel_surf.^2);
fVel_Elmer (mask==0) = NaN;
fVel_Elmer_diff = abs(fVel_Elmer - velo2011)./abs(velo2011)*100;
%fVel_Elmer_diff = fVel_Elmer - velo2011;
fVel_Elmer_diff(mask==0) = nan;
%fVel_Elmer_diff(velo2011==0) = nan;
%fVel_ini_Elmer = sqrt(Elmer_xvel_ini_surf.^2 + Elmer_yvel_ini_surf.^2+Elmer_zvel_ini_surf.^2);
%fVel_ini_Elmer (mask==0) = NaN;
%Elmer_fVel_change = fVel_Elmer - fVel_ini_Elmer;

beta_Elmer = log10(Elmer_beta_surf);
beta_Elmer (mask==0) = NaN;
% beta_ini_Elmer = log10(Elmer_beta_ini_surf);
% beta_ini_Elmer (mask==0) = NaN;
%*************************************************************************%

fprintf('Plotting... \n');

figure;
H=imagesc(fVel_Elmer_diff);colorbar;colormap(jet);%title('Elmer: Vel_{surf}-Vel_{surf-control}');
hold on;
%set(H,'alphadata',mask~=0);
%set(H,'alphadata',~isnan(fVel_Elmer_diff));
set(H,'alphadata',~isnan(velo2011));
    %plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    %plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
    plot(B3_side_FoPr(1,:),B3_side_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
    %plot(I_xllarge:435,linspace(I_yularge,I_yularge,size(I_xllarge:435,2)),'-.m','LineWidth',2),hold on;
    %plot(linspace(I_xllarge,I_xllarge,size(I_yularge:309,2)),I_yularge:309,'-.m','LineWidth',2),hold on;
    %plot(I_xl:427,linspace(I_yu,I_yu,size(I_xl:427,2)),'-.y','LineWidth',2),hold on;
    %plot(linspace(I_xl,I_xl,size(I_yu:310,2)),I_yu:310,'-.y','LineWidth',2),hold on;
%     plot(proN.ArrayX(1).data,proN.ArrayY(1).data,'--w','LineWidth',2);
%     plot(proM.ArrayX(1).data,proM.ArrayY(1).data,'--w','LineWidth',2);
%     plot(proS.ArrayX(1).data,proS.ArrayY(1).data,'--w','LineWidth',2);
    edittingASFplot(xaxi,yaxi,colorbarlim(count,:),xlab,ylab,0,0,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
    %colormap(darkb2r(colorbarlim(count,1),colorbarlim(count,2)));
%saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_fvel.tif']),'tif');
count = count + 1;

end
end

fprintf('The end \n');
 


