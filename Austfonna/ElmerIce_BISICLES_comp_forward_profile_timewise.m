close all;
clear all;
fprintf('Loading data... \n');
path2 = 'J:\ASF_data\B_results\';
load ASF_B3_mask.mat;
topg = ncread([path2,'ASF95Efor_BISICLES.nc'],'topg');
topg = rot90(topg);
thk = ncread([path2,'ASF95Efor_BISICLES.nc'],'thk');
thk = rot90(thk);
%*************************************************************************%

fprintf('Preparing data... \n');
obvel2011x = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'xvel');
obvel2011x = rot90(obvel2011x);
obvel2011y = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'yvel');
obvel2011y = rot90(obvel2011y);
velo2011 = sqrt(obvel2011x.^2 + obvel2011y.^2);
velo2011(mask == 0)= NaN;

base_name_Elmer = 'Elmer';
base_name_BISICLES = 'BISICLES';
time_name = {'95to11'};
smb_name = {'elevcorr_ini','elevcorr_yr05','elevcorr_yr10','elevcorr'};
beta_name = {'celmer_linear'};
beta_legend_name = {'1995Jan','2005Jan','2010Jan','2011Dec'};
time_legend_name = {'','2005','2010','2011'};
temp_name = {'temp95'};
E_beta_ini_name = {'nothermal_ebet'};
E_path_name = 'J:\ASF_data\E_results';
B_path_name = 'J:\ASF_data\B_results';

%*****************************configure***********************************%
  grey=[0.8,0.8,0.8]; 
  %linecolor= ['c', 'm','b','g'];
  x_B3 = [322 352 377 413 419 426 430];
  y_B3 = [188 248 265 270 268 265 259];
  x_lab_B3 = ((340:20:420)-322)*0.4;
  xlimit_B3 = [322 430];
  yy1limit_B3 = [-130 1000]; 
  %datalimit_B3 = [0 800]; 
  %datalimit_B3 = [-200 800]; 
  %betalimit_B3 = [-9 0];
  yy2limit_B3 = [-50 40];
  xrange = 'x_B3'; yrange = 'y_B3';
  xlimit = 'xlimit_B3'; yy1limit = 'yy1limit_B3';
  x_lab = 'x_lab_B3';yy2limit = 'yy2limit_B3';
  color = {'r','b'};
  linestyle = {'-','--','-',':'};
  %legend_posi =  [0.51,0.6,0.1,0.3];
  legend_posi = 'Best';
  counter = 1;
  [ArrayXG(1).data,ArrayYG(1).data,ArrayCG(1).data]=improfile(topg,eval(xrange),...
     eval(yrange),'bicubic');
  [ArrayXG(2).data,ArrayYG(2).data,ArrayCG(2).data]=improfile(topg+thk,eval(xrange),...
     eval(yrange),'bicubic');
  
  
 figure(1);
  fillX = [ArrayXG(1).data; flipud(ArrayXG(1).data)];
  fillY = [ArrayCG(2).data; flipud(ArrayCG(1).data)];
  Hfill = fill(fillX,fillY,grey);hold on;
  set(Hfill,'LineWidth',1);   
  xlabel('Distance(km)','FontSize',15);
  
%*****************************BISICLES************************************%
for j = 1: length(beta_name)

for i = 1: length(smb_name)
%*******************************Elmer**************************************%

BISICLES_xbVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.nc']),'xbVel'));
BISICLES_ybVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.nc']),'ybVel'));
BISICLES_xfVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.nc']),'xfVel'));
BISICLES_yfVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.nc']),'yfVel'));
BISICLES_FS= rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.nc']),'Z_surface'));
BISICLES_beta = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.nc']),'basal_friction'));
beta_BISICLES =BISICLES_beta.*(10^-6);
beta_BISICLES =log10(BISICLES_beta.*(10^-6));
beta_BISICLES(mask==0) = NaN;

BISICLES_xbVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_elevcorr_ini.nc']),'xbVel'));
BISICLES_ybVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_elevcorr_ini.nc']),'ybVel'));
BISICLES_xfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_elevcorr_ini.nc']),'xfVel'));
BISICLES_yfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_elevcorr_ini.nc']),'yfVel'));
BISICLES_FS_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_elevcorr_ini.nc']),'Z_surface'));

% BISICLES_xbVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'xbVel'));
% BISICLES_ybVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'ybVel'));
% BISICLES_xfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'xfVel'));
% BISICLES_yfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'yfVel'));
%BISICLES_FS_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_20thsmb.nc']),'Z_surface'));
fVel_BISICLES = sqrt(BISICLES_xfVel.^2 + BISICLES_yfVel.^2);
fVel_BISICLES (mask==0) = NaN;
fVel_ini_BISICLES = sqrt(BISICLES_xfVel_ini.^2 + BISICLES_yfVel_ini.^2);
fVel_ini_BISICLES (mask==0) = NaN;
BISICLES_fVel_change = fVel_BISICLES ;%- fVel_ini_BISICLES;

BISICLES_FS_ini (mask==0) = NaN;
BISICLES_FS (mask==0) = NaN;
BISICLES_FS_change = BISICLES_FS - BISICLES_FS_ini;

%*******************************Elmer**************************************%

filename = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_',cell2mat(smb_name(i)),'.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,cell2mat(E_beta_ini_name(i)),'.ep.mat'];
filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(j)),'_elevcorr_ini.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_','20thsmb.ep9.mat'];
load(cell2mat(filename));
Eini=load(cell2mat(filename_ini));

 
Elmer_xvel_bed = matrixdata(2).data_layer(1).data;
Elmer_yvel_bed = matrixdata(3).data_layer(1).data;
Elmer_zvel_bed = matrixdata(4).data_layer(1).data;
Elmer_xvel_surf = matrixdata(2).data_layer(2).data;
Elmer_yvel_surf = matrixdata(3).data_layer(2).data;
Elmer_zvel_surf = matrixdata(4).data_layer(2).data;
Elmer_FS_surf = matrixdata(9).data_layer(1).data;
data_Ebeta = matrixdata(8).data_layer(1).data;
Elmer_beta = data_Ebeta;
data_Ebeta = log10(data_Ebeta);
% Elmer_xvel_ini_bed = Eini.matrixdata(2).data_layer(1).data;
% Elmer_yvel_ini_bed = Eini.matrixdata(3).data_layer(1).data;
% Elmer_zvel_ini_bed = Eini.matrixdata(4).data_layer(1).data;
Elmer_xvel_ini_surf = Eini.matrixdata(2).data_layer(2).data;
Elmer_yvel_ini_surf = Eini.matrixdata(3).data_layer(2).data;
Elmer_zvel_ini_surf = Eini.matrixdata(4).data_layer(2).data;
Elmer_FS_ini_surf = Eini.matrixdata(9).data_layer(1).data;
% 
 Elmer_FS_ini_surf (mask==0) = NaN;
Elmer_FS_surf (mask==0) = NaN;
Elmer_FS_change = Elmer_FS_surf-Elmer_FS_ini_surf;

bVel_Elmer = sqrt(Elmer_xvel_bed.^2 +Elmer_yvel_bed.^2+Elmer_zvel_bed.^2);
bVel_Elmer (mask==0) = NaN;

% bVel_ini_Elmer = sqrt(Elmer_xvel_ini_bed.^2 + Elmer_yvel_ini_bed.^2+Elmer_zvel_ini_bed.^2);
% bVel_ini_Elmer (mask==0) = NaN;
% Elmer_bVel_change = bVel_Elmer - bVel_ini_Elmer;

fVel_Elmer = sqrt(Elmer_xvel_surf.^2 +Elmer_yvel_surf.^2+Elmer_zvel_surf.^2);
fVel_Elmer (mask==0) = NaN;

fVel_ini_Elmer = sqrt(Elmer_xvel_ini_surf.^2 + Elmer_yvel_ini_surf.^2+Elmer_zvel_ini_surf.^2);
fVel_ini_Elmer (mask==0) = NaN;
Elmer_fVel_change = fVel_Elmer;% - fVel_ini_Elmer;
 %   *************************************************************************
  fprintf('plotting profile...\n');
  % elevation change
 [ArrayX(1).data,ArrayY(1).data,ArrayC(1).data]=improfile(Elmer_FS_change,eval(xrange),...
     eval(yrange),'bicubic');
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data]=improfile(BISICLES_FS_change,eval(xrange),...
     eval(yrange),'bicubic');

[ArrayX(3).data,ArrayY(3).data,ArrayC(3).data]=improfile(Elmer_fVel_change,eval(xrange),...
    eval(yrange),'bicubic');

 [ArrayX(4).data,ArrayY(4).data,ArrayC(4).data]=improfile(BISICLES_fVel_change,eval(xrange),...
     eval(yrange),'bicubic');
  for k =1:2
      if i==1
          LC = 'k';
      else
          LC= cell2mat(color(k));
      end
    
    figure(1);
    H(counter) = plot(ArrayX(k+2).data,ArrayC(k+2).data,...
         cell2mat(linestyle(i)),'LineWidth',2,'Color',LC);
    hold on;
    
    %beta_legend_name_input = cell2mat(beta_legend_name(j));
    time_legend_name_input = cell2mat(time_legend_name(i));
  
    %legend_vel(counter) = {[time_legend_name_input,': |u|_{',beta_legend_name_input,'}']};
    switch k
        case 1
        %legend_FS(counter) = {'Elmer/Ice:\DeltaFS'};
        legend_vel(counter) = {[time_legend_name_input,': |u|_{Elmer/Ice}']};
        case 2
        %legend_FS(counter) = {'BISICLES:\DeltaFS'};
        legend_vel(counter) = {[time_legend_name_input,': |u|_{BISICLES}']};
    end
    counter = counter + 1;
  end
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','c','LineStyle','-.','LineWidth',1.5);hold on;
  
end
for k =1:2
 if 0  
  if i==1
  else
    LC= cell2mat(color(k));
    figure(2),
     plot(ArrayX(k).data,ArrayC(k).data,...
         cell2mat(linestyle(i)),'LineWidth',2,'Color',LC);
    hold on;
    
   
  end
end
end
end

set(gca,'XTickLabel',eval(x_lab),...
    'XLim',eval(xlimit), 'YLim',eval(yy1limit),...
    'FontSize',15);%somehow the y direction of grid has been reversed
    ylabel('Elevtion (m) or Speed (m a^{-1})');
    ylim(eval(yy1limit));
% axis tight;
    % 'YTick',linspace(yy1limit_B3(1),yy1limit_B3(2),11), ...
if 0
   LH=legend(gca(figure(1)),[Hfill;H(1);H(3);H(4);H(5);H(6);H(7);H(8);lineH], ...
       'Glacier_{initial}',...
       '|u|_{initial}',...
       cell2mat(legend_vel(3)),cell2mat(legend_vel(4)),cell2mat(legend_vel(5)),...
       cell2mat(legend_vel(6)),cell2mat(legend_vel(7)),cell2mat(legend_vel(8)),...
       'sea level',...
       'Location',legend_posi);
    set(LH,'color','w','FontSize',30);
end
if 0
set(gca(figure(2)),...
    'XTickLabel',eval(x_lab),...
    'XLim',eval(xlimit), 'YLim',eval(yy2limit),...
    'FontSize',15);%somehow the y direction of grid has been reversed
    ylabel('Elevtion Change (m)');xlabel('Distance(km)','FontSize',15);
    ylim(eval(yy2limit));      
    %axis tight;
end
save('profile_B3.mat','ArrayX','ArrayY');


