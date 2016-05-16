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
  yy2limit_B3 = [-100 40];
  xrange = 'x_B3'; yrange = 'y_B3';
  xlimit = 'xlimit_B3'; yy1limit = 'yy1limit_B3';
  x_lab = 'x_lab_B3';yy2limit = 'yy2limit_B3';
  smb_color = {'y','r','c'};
  %legend_posi =  [0.51,0.6,0.1,0.3];
  legend_posi = 'NorthEastOutside';
  counter = 1;
  [ArrayXG(1).data,ArrayYG(1).data,ArrayCG(1).data]=improfile(topg,eval(xrange),...
     eval(yrange),'bicubic');
  [ArrayXG(2).data,ArrayYG(2).data,ArrayCG(2).data]=improfile(topg+thk,eval(xrange),...
     eval(yrange),'bicubic');
  
  
%*****************************BISICLES************************************%
for j = 1: length(smb_name)
  figure,%subplot(2,1,i)
  fillX = [ArrayXG(1).data; flipud(ArrayXG(1).data)];
  fillY = [ArrayCG(2).data; flipud(ArrayCG(1).data)];
  Hfill(1) = fill(fillX,fillY,grey);hold on;
  set(Hfill(1),'LineWidth',2);
  xlabel('Distance(km)','FontSize',15);
for i = 1: length(beta_name)
    
BISICLES_xbVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'xbVel'));
BISICLES_ybVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'ybVel'));
BISICLES_xfVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'xfVel'));
BISICLES_yfVel = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'yfVel'));
BISICLES_FS= rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'Z_surface'));
BISICLES_beta = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.nc']),'basal_friction'));
beta_BISICLES =BISICLES_beta.*(10^-6);
beta_BISICLES =log10(BISICLES_beta.*(10^-6));
beta_BISICLES(mask==0) = NaN;

BISICLES_xbVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'xbVel'));
BISICLES_ybVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'ybVel'));
BISICLES_xfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'xfVel'));
BISICLES_yfVel_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'yfVel'));
BISICLES_FS_ini = rot90(ncread(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.nc']),'Z_surface'));

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

filename = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,cell2mat(E_beta_ini_name(i)),'.ep.mat'];
filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_elevcorr_ini.ep10.mat'];
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
  
 [ArrayX(1).data,ArrayY(1).data,ArrayC(1).data]=improfile(Elmer_FS_change,eval(xrange),...
     eval(yrange),'bicubic');
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data]=improfile(BISICLES_FS_change,eval(xrange),...
     eval(yrange),'bicubic');
%  [ArrayX(3).data,ArrayY(3).data,ArrayC(3).data]=improfile(fVel_Elmer,eval(xrange),...
%      eval(yrange),'bicubic');
[ArrayX(3).data,ArrayY(3).data,ArrayC(3).data]=improfile(Elmer_fVel_change,eval(xrange),...
    eval(yrange),'bicubic');
%  [ArrayX(3).data,ArrayY(3).data,ArrayC(3).data]=improfile(data_Ebeta,eval(xrange),...
%      eval(yrange),'bicubic');
%  [ArrayX(4).data,ArrayY(4).data,ArrayC(4).data]=improfile(fVel_BISICLES,eval(xrange),...
%      eval(yrange),'bicubic');
 [ArrayX(4).data,ArrayY(4).data,ArrayC(4).data]=improfile(BISICLES_fVel_change,eval(xrange),...
     eval(yrange),'bicubic');
%   Hfill(counter + 1) = plot(ArrayX(1).data,ArrayC(1).data);hold on;
%   set(Hfill(counter + 1),'LineWidth',2);
%   xlabel('Distance along section (km)','FontSize',15);
  
  for k=1:2   
    [AX,Hplot1(counter).data,Hplot2(counter).data] = plotyy(ArrayX(k+2).data,ArrayC(k+2).data,ArrayX(k).data,ArrayC(k).data,'plot');
    hold on;
    beta_legend_name_input = cell2mat(beta_legend_name(j));
    
    %legend_beta(counter) = {[' \beta_{',beta_legend_name_input,'}']};
    
    
    set(AX(1),'YLim',eval(yy1limit),'XLim',eval(xlimit),'YColor','k',...
     'XTickLabel',eval(x_lab),'FontSize',15,'ytick',linspace(yy1limit_B3(1),yy1limit_B3(2),11));
    set(get(AX(1),'Ylabel'),'string','Elevtion (m) Speed (m a^{-1})','FontSize',15); 
    %ylabel('Elevtion (m), or Speed (m/a)');
    set(AX(2),'YLim',eval(yy2limit),'XLim',eval(xlimit),'YColor','k',...
      'XTickLabel',eval(x_lab),'FontSize',15,'ytick',linspace(-1*yy2limit_B3(2),yy2limit_B3(2),9));
    set(get(AX(2),'Ylabel'),'string','Elevtion Change (m)','FontSize',15); 
%     if i==3||i==4
%         
%         set(Hplot2(counter).data,'Color','none','LineStyle', '--','LineWidth',1.5);
%     else   
    switch k
        case 1
        set(Hplot1(counter).data,'LineStyle', '-','LineWidth',2,'Color','c');  
        set(Hplot2(counter).data,'Color',get(Hplot1(counter).data,'Color'),'LineStyle', '--','LineWidth',2);
        legend_FS(counter) = {'Elmer/Ice:\DeltaFS'};
        legend_vel(counter) = {'Elmer/Ice:|u|'};
        case 2
        set(Hplot1(counter).data,'LineStyle', '-','LineWidth',2,'Color','r');  
        set(Hplot2(counter).data,'Color',get(Hplot1(counter).data,'Color'),'LineStyle', '--','LineWidth',2);
        legend_FS(counter) = {'BISICLES:\DeltaFS'};
        legend_vel(counter) = {'BISICLES:|u|'};
    end
    counter = counter + 1;
  end
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','b','LineStyle','-.','LineWidth',1.5);hold on;
     %ylabel('\beta (MPa/a/m)');
     %set(gca,'XTickLabel',eval(x_lab),'FontSize',15);
     %figure,imagesc(data_Ebeta),colorbar,colormap(jet);
end
end

LH=legend(gca(figure(1)),[Hfill;...
       Hplot1(7).data;Hplot1(8).data;...
       Hplot2(7).data;Hplot2(8).data;...
       lineH],...
       'Glacier_{initial}',...
       cell2mat(legend_vel(7)),cell2mat(legend_vel(8)),...
       cell2mat(legend_FS(7)),cell2mat(legend_FS(8)),...
       'sea level',...
       'Location',legend_posi);
    set(LH,'color','w','FontSize',30);
if 0
   LH=legend([Hfill;...
       Hplot1(1).data;Hplot1(2).data;Hplot1(3).data;Hplot1(4).data;...
       Hplot1(5).data;Hplot1(6).data;Hplot1(7).data;Hplot1(8).data;...
       Hplot2(1).data;Hplot2(2).data;Hplot2(3).data;Hplot2(4).data;...
       Hplot2(5).data;Hplot2(6).data;Hplot2(7).data;Hplot2(8).data;...
       lineH],...
       'Glacier_{initial}',...
       cell2mat(legend_vel(1)),cell2mat(legend_vel(2)),cell2mat(legend_vel(3)),cell2mat(legend_vel(4)),...
       cell2mat(legend_vel(5)),cell2mat(legend_vel(6)),cell2mat(legend_vel(7)),cell2mat(legend_vel(8)),...
       cell2mat(legend_FS(1)),cell2mat(legend_FS(2)),cell2mat(legend_FS(3)),cell2mat(legend_FS(4)),...
       cell2mat(legend_FS(5)),cell2mat(legend_FS(6)),cell2mat(legend_FS(7)),cell2mat(legend_FS(8)),...
       'sea level',...
       'Location',legend_posi);
    set(LH,'color','w','FontSize',15);
end
       %      'Glacier_{initial}','FS_E_l_m_e_r_/_I_c_e','FS_B_I_S_I_C_L_E_S',...
%      '\Deltau_E_l_m_e_r_/_I_c_e','lg(\beta_E_l_m_e_r_/_I_c_e)','\Deltau_B_I_S_I_C_L_E_S',' ','water line',...
%      'Location',legend_posi);
% set(LH,'color','none')
% data_B = BISICLES_beta(mask==2);
% data_E = Elmer_beta(mask==2);
% % data_B = fVel_BISICLES;
% % data_E = fVel_Elmer;
% data_B = data_B.*(10^-6);
% data_B(isnan(data_B))=0;
% data_E(isnan(data_E))=0;
% ErrC =  (sum(abs(data_B-data_E))/sum(abs(data_E)))*100

save('profile_B3.mat','ArrayX','ArrayY');



