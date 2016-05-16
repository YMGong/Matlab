
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


%*************************************************************************%

fprintf('Preparing data... \n');
obvel2011x = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'xvel');
obvel2011x = rot90(obvel2011x);
obvel2011y = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'yvel');
obvel2011y = rot90(obvel2011y);
velo2011 = sqrt(obvel2011x.^2 + obvel2011y.^2);
velo2011(mask == 0)= NaN;
%velo2011(velo2011<3) = NaN;

base_name_Elmer = 'Elmer';
base_name_BISICLES = 'BISICLES';
%volume_name = {'asf','b3'};
time_name = {'95to11'};
%smb_name = {'elevcorr','noelevcorr','20thsmb'};
smb_name = {'elevcorr'};
%E_smb_name = {'elecorr'};
%beta_name = {'celmer_linear','celmer_stepc95','celmer_stepc2011'};
beta_name = {'celmer_linear'};
%beta_name = {'linearbeta_','95beta_','2011beta_'};
temp_name = {'temp95'};
B_beta_ini_name = {'20yr'};
E_beta_ini_name = {'nothermal_ebet'};
B_path_name = '../B_results';
E_path_name = '../E_results';

%*****************************configure***********************************%
%xaxi = [150 500];yaxi = [50 370];
xaxi = [305 445];yaxi = [180 320];
xlab = (15:20:135)*0.400; ylab = (140:-20:0)*0.400;  
x_tick = 0;
y_tick = 0;

gridon = 1; lineon = 0;boxon = 1;
axi_units = 'km';
grey = [0.7,0.7,0.7];
%*****************************BISICLES************************************%

for i = 1: length(beta_name)
for j = 1: length(smb_name)
BISICLES_xbVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'xbVel'));
BISICLES_ybVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'ybVel'));
BISICLES_xfVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'xfVel'));
BISICLES_yfVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'yfVel'));
BISICLES_thickness= rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'thickness'));
BISICLES_beta = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'basal_friction'));

BISICLES_xbVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'xbVel'));
BISICLES_ybVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'ybVel'));
BISICLES_xfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'xfVel'));
BISICLES_yfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'yfVel'));
BISICLES_thickness_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'thickness'));
% BISICLES_xbVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp2011_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'xbVel'));
% BISICLES_ybVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp2011_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'ybVel'));
% BISICLES_xfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp2011_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'xfVel'));
% BISICLES_yfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp2011_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'yfVel'));
% BISICLES_thickness_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp2011_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'thickness'));
% BISICLES_xbVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'xbVel'));
% BISICLES_ybVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'ybVel'));
% BISICLES_xfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'xfVel'));
% BISICLES_yfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'yfVel'));
% BISICLES_thickness_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'thickness'));
% BISICLES_xbVel_ini = rot90(ncread([B_path_name,'/',base_name_BISICLES,cell2mat(B_beta_ini_name(i)),'.nc'],'xbVel'));
% BISICLES_ybVel_ini= rot90(ncread([B_path_name,'/',base_name_BISICLES,cell2mat(B_beta_ini_name(i)),'.nc'],'ybVel'));
% BISICLES_xfVel_ini = rot90(ncread([B_path_name,'/',base_name_BISICLES,cell2mat(B_beta_ini_name(i)),'.nc'],'xfVel'));
% BISICLES_yfVel_ini= rot90(ncread([B_path_name,'/',base_name_BISICLES,cell2mat(B_beta_ini_name(i)),'.nc'],'yfVel'));
% BISICLES_thickness_ini= rot90(ncread([B_path_name,'/',base_name_BISICLES,cell2mat(B_beta_ini_name(i)),'.nc'],'thickness'));

BISICLES_thickness_ini (mask==0) = NaN;
BISICLES_thickness (mask==0) = NaN;
BISICLES_thickness_change = BISICLES_thickness - BISICLES_thickness_ini;

% 
% Vel_BISICLES = sqrt(BISICLES_xVel.^2 + BISICLES_yVel.^2+BISICLES_zVel.^2);
% Vel_BISICLES (mask==0) = NaN;
% Vel_ini_BISICLES = sqrt(BISICLES_xVel_ini.^2 + BISICLES_yVel_ini.^2+BISICLES_zVel_ini.^2);
% Vel_ini_BISICLES (mask==0) = NaN;
% BISICLES_Vel_change = Vel_BISICLES - Vel_ini_BISICLES;


bVel_BISICLES = sqrt(BISICLES_xbVel.^2 + BISICLES_ybVel.^2);
bVel_BISICLES (mask==0) = NaN;

% bVel_ini_BISICLES = sqrt(BISICLES_xbVel_ini.^2 + BISICLES_ybVel_ini.^2);
% bVel_ini_BISICLES (mask==0) = NaN;
% BISICLES_bVel_change = bVel_BISICLES - bVel_ini_BISICLES;

fVel_BISICLES = sqrt(BISICLES_xfVel.^2 + BISICLES_yfVel.^2);
fVel_BISICLES (mask==0) = NaN;

% fVel_ini_BISICLES = sqrt(BISICLES_xfVel_ini.^2 + BISICLES_yfVel_ini.^2);
% fVel_ini_BISICLES (mask==0) = NaN;
%BISICLES_fVel_change = fVel_BISICLES - fVel_ini_BISICLES;

%beta_BISICLES =log10(BISICLES_beta.*(10^-6));
beta_BISICLES =BISICLES_beta.*(10^-6);
beta_BISICLES(mask==0) = NaN;

%*******************************Elmer**************************************%

filename = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,cell2mat(E_beta_ini_name(i)),'.ep.mat'];
% filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_','20thsmb.ep10.mat'];
filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_','ini.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_temp2011_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.ep10.mat'];

load(cell2mat(filename));
Eini=load(cell2mat(filename_ini));

Elmer_xvel_bed = matrixdata(2).data_layer(1).data;
Elmer_yvel_bed = matrixdata(3).data_layer(1).data;
Elmer_zvel_bed = matrixdata(4).data_layer(1).data;
Elmer_xvel_surf = matrixdata(2).data_layer(2).data;
Elmer_yvel_surf = matrixdata(3).data_layer(2).data;
Elmer_zvel_surf = matrixdata(4).data_layer(2).data;
Elmer_height_surf = matrixdata(1).data_layer(1).data;
Elmer_beta_surf = matrixdata(8).data_layer(1).data;

% Elmer_xvel_ini_bed = Eini.matrixdata(2).data_layer(1).data;
% Elmer_yvel_ini_bed = Eini.matrixdata(3).data_layer(1).data;
% Elmer_zvel_ini_bed = Eini.matrixdata(4).data_layer(1).data;
% Elmer_xvel_ini_surf = Eini.matrixdata(2).data_layer(2).data;
% Elmer_yvel_ini_surf = Eini.matrixdata(3).data_layer(2).data;
% Elmer_zvel_ini_surf = Eini.matrixdata(4).data_layer(2).data;
 Elmer_height_ini_surf = Eini.matrixdata(1).data_layer(1).data;

Elmer_height_ini_surf (mask==0) = NaN;
Elmer_height_surf (mask==0) = NaN;
Elmer_thickness_change = Elmer_height_surf - Elmer_height_ini_surf;

bVel_Elmer = sqrt(Elmer_xvel_bed.^2 +Elmer_yvel_bed.^2+Elmer_zvel_bed.^2);
bVel_Elmer (mask==0) = NaN;

% bVel_ini_Elmer = sqrt(Elmer_xvel_ini_bed.^2 + Elmer_yvel_ini_bed.^2+Elmer_zvel_ini_bed.^2);
% bVel_ini_Elmer (mask==0) = NaN;
% Elmer_bVel_change = bVel_Elmer - bVel_ini_Elmer;

fVel_Elmer = sqrt(Elmer_xvel_surf.^2 +Elmer_yvel_surf.^2+Elmer_zvel_surf.^2);
fVel_Elmer (mask==0) = NaN;

% fVel_ini_Elmer = sqrt(Elmer_xvel_ini_surf.^2 + Elmer_yvel_ini_surf.^2+Elmer_zvel_ini_surf.^2);
% fVel_ini_Elmer (mask==0) = NaN;
% Elmer_fVel_change = fVel_Elmer - fVel_ini_Elmer;

beta_Elmer = log10(Elmer_beta_surf);
beta_Elmer (mask==0) = NaN;

%*************************************************************************%

fprintf('Plotting... \n');
% thickness don't need it now
if 0
figure;
colorbarlim=[0,538]; colorbar_lab = (0:1)*100;c_units = 'm';

H=imagesc(BISICLES_thickness);hold on;colorbar;
title('thickness-BISICLES');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_thickness.tif']),'tif');


figure;
H=imagesc(Elmer_height_surf);hold on;colorbar;title('thickness-Elmer');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
 
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_thickness.tif']),'tif');

% figure;
% imagesc((thickness-Elmer_height_surf));colorbar;title('thickness:BISICLES - Elmer');
% colormap(darkb2r(-264,130));
% saveas(gcf,cell2mat(['.\pics\Diff_',cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_thickness.tif']),'tif');

colorbarlim=[-20,20]; colorbar_lab = (0:1)*100;c_units = 'm';
figure;
H=imagesc(BISICLES_thickness_change);hold on;colorbar;title('thickness-change-BISICLES');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');

if strcmp(cell2mat(beta_name(i)),'95beta_')
    colormap(darkb2r(-5,20));
else
    colormap(darkb2r(-15,15));
end

saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_thickness_change.tif']),'tif');

colorbarlim=[-20,20];
figure;
H=imagesc(Elmer_thickness_change);hold on;colorbar;title('thickness-change-Elmer');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
colormap(darkb2r(-50,50));
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_thickness_change.tif']),'tif');
%*************************************************************************%

% basal velocity do not need now

figure;
colorbarlim=[0,700]; colorbar_lab = (0:1)*100;c_units = 'm/yr';

H=imagesc(bVel_BISICLES);hold on;colorbar;title('bvel-BISICLES');
%set(gca,'Clim',[0.0,600]);  
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_bvel.tif']),'tif');

figure;
H=imagesc(bVel_Elmer);hold on;colorbar;title('bvel-Elmer');
%set(gca,'Clim',[0.0,600]);  
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_bvel.tif']),'tif');


colorbarlim=[-100,50];
figure;
H=imagesc((bVel_BISICLES-bVel_Elmer));hold on;colorbar;title('bvel:BISICLES - Elmer');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
colormap(darkb2r(-200,80));
saveas(gcf,cell2mat(['.\pics\Diff_',cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_bvel.tif']),'tif');
end
colorbarlim=[-3,3];colorbar_lab = (0:1)*100;c_units = 'm/yr';
if 0
figure;
H=imagesc(BISICLES_bVel_change);hold on;colorbar;title('bvel-change-BISICLES');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
%colormap(darkb2r(-100,400));
saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_bvel_change.tif']),'tif');

figure;
H=imagesc(Elmer_bVel_change);hold on;colorbar;title('bvel-change-Elmer');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
%colormap(darkb2r(-100,400));
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_bvel_change.tif']),'tif');
end
%*************************************************************************%
% surface velocity only need the absolute
colorbar_lab = (0:1)*100;c_units = 'm a^{-1}';
colorbarlim=[0,350];
if 0
figure;
H=imagesc(fVel_BISICLES);hold on;colorbar;title('fVel-BISICLES');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_fvel.tif']),'tif');

figure;
H=imagesc(fVel_Elmer);colorbar;title('fVel-Elmer');hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_fvel.tif']),'tif');


colorbarlim=[-3,3];
figure;
H=imagesc((fVel_BISICLES-fVel_Elmer));colorbar;title('fvel:BISICLES - Elmer');hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
colormap(darkb2r(-30,30));
saveas(gcf,cell2mat(['.\pics\Diff_',cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_fvel.tif']),'tif');
end

colorbarlim=[0,1000];
%colorbarlim=[-30,90];

figure;
H=imagesc(fVel_BISICLES);colorbar;%title('BISICLES: Vel_{surf}-Vel_{surf-control}');
hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
    plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
%colormap(darkb2r(-30,200));
saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_fvel_change.tif']),'tif');
figure;
H=imagesc(fVel_Elmer);colorbar;%title('Elmer: Vel_{surf}-Vel_{surf-control}');
hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
    plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
%colormap(darkb2r(0,500));
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_fvel_change.tif']),'tif');
figure;
H=imagesc(fVel_BISICLES - fVel_Elmer);colorbar;%title('Elmer: Vel_{surf}-Vel_{surf-control}');
hold on;
colorbarlim=[-100,100];gridon = 0; lineon = 0;boxon = 0;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
    plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
%axes('xcolor','w','ycolor','w');
colormap(darkb2r(-100,100));
%*************************************************************************%
% beta do not need now
if 0
colorbar_lab = (0:1)*100;c_units = 'MPa/a/m';
colorbarlim=[-6,0];
figure;
H=imagesc(beta_BISICLES);hold on;colorbar;title('beta-BISICLES');
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_BISICLES,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_beta.tif']),'tif');

figure;
H=imagesc(beta_Elmer);colorbar;title('beta-Elmer');hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
saveas(gcf,cell2mat(['.\pics\',base_name_Elmer,cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_beta.tif']),'tif');


colorbarlim=[-1,1];
figure;
H=imagesc((beta_BISICLES-beta_Elmer));colorbar;title('beta:BISICLES - Elmer');hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
colormap(darkb2r(-1,1));
saveas(gcf,cell2mat(['.\pics\Diff_',cell2mat(beta_name(i)),time_name(1),'_',cell2mat(smb_name(j)),'_beta.tif']),'tif');
end
% figure(16);
% imagesc(Vel_BISICLES);colorbar;title('Vel-BISICLES');
% set(gca,'Clim',[0.0,800]);  
% 
% figure(17);
% imagesc(BISICLES_Vel_change);colorbar;title('Vel-change-BISICLES');
% colormap(darkb2r(-300,500));

end
if 0
%checking B_C do not need now
%currentT = 20.0021;
%currentT =  10.0061;
currentT = 5.0081;
Tini = 0.0101;
Tinterv = 21.0;
colorbar_lab = (0:1)*100;c_units = 'MPa/a/m';
beta_final=ncread('/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/frc.velomatch.monthly000022.new.2d.nc','C');
beta_final=rot90(beta_final.*(10^-6));
beta_final(mask==0) = NaN;
beta_initial=ncread('/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/frc.velomatch.monthly000000.new.2d.nc','C');
beta_initial=rot90(beta_initial.*(10^-6));
beta_initial(mask==0) = NaN;

beta_tmp = (1 - (floor(currentT - Tini) / Tinterv)) * beta_initial +...
    (floor(currentT - Tini) / Tinterv) * beta_final;

figure;
H=imagesc((beta_BISICLES-beta_tmp));colorbar;title('beta:BISICLES - Elmer');hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,[-1 1],xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
colormap(darkb2r(-1,1));


currentTD = ((beta_BISICLES-beta_initial)./(beta_final-beta_initial)).*Tinterv;
slope = (beta_final-beta_initial)./Tinterv;
currentSlope = (beta_BISICLES-beta_initial)./(currentT-Tini);

% Tintern = floor(currentT - Tini) ./ data;
% Tdiff = data.*Tinterv;

figure,imagesc(slope),colorbar,colormap(jet);title('slope');
 edittingASFplot(xaxi,yaxi,[-0.1 0.1],xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
figure,imagesc(currentTD),colorbar,colormap(jet);title('currentTD');
figure,imagesc(currentSlope),colorbar,colormap(jet);title('currentSlope');
 edittingASFplot(xaxi,yaxi,[-0.1 0.1],xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');

figure;
H=imagesc((fVel_Elmer-fVel_BISICLES));colorbar;title('velocity:2011 observation - BISICLES');hold on;
set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    edittingASFplot(xaxi,yaxi,[-80 80],xlab,ylab,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
colormap(darkb2r(-80,80));
end

end
fprintf('The end \n');
 


