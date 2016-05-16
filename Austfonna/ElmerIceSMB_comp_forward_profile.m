close all;
clear all;
fprintf('Loading data... \n');

path2 = 'J:\ASF_data\B_results\';
addpath 'J:\ASF_data\hirham';
load ASF_B3_mask.mat;
smbini=load('../E_results/Elevcorr_smb_ini.mat');
smbfin=load('../E_results/Elevcorr_smb_fin.mat');
topg = ncread([path2,'ASF95Efor_BISICLES.nc'],'topg');
topg = rot90(topg);
thk = ncread([path2,'ASF95Efor_BISICLES.nc'],'thk');
thk = rot90(thk);
grey=[0.6,0.6,0.6]; 

% smbini_data = smbini.matrixdata(1).data_layer.data/910.0;
% smbini_data(mask==0)=nan;
% smbfin_data = smbfin.matrixdata(1).data_layer.data/910.0;
% smbfin_data(mask==0)=nan;
smbini=0.0;
for i = 60:71
smbini = smbini + ncread(['hirham.smb.jan.1990.jan.2011.elevcorrec0000',num2str(i),'.2d.nc'],'smb');
end
smbini = rot90(smbini);
smbini(mask==0)=nan;

smbfin=0.0;
for i = 252:263
smbfin = smbfin + ncread(['hirham.smb.jan.1990.jan.2011.elevcorrec000',num2str(i),'.2d.nc'],'smb');
end
smbfin = rot90(smbfin);
smbfin(mask==0)=nan;

% figure,imagesc(smbini),colormap(jet),colorbar;
% figure,imagesc(smbfin),colormap(jet),colorbar;

smbininocorr=0.0;
for i = 60:71
smbininocorr = smbininocorr + ncread(['hirham.smb.jan.1990.jan.2011.nocorrec0000',num2str(i),'.2d.nc'],'smb');
end
smbininocorr = rot90(smbininocorr);
smbininocorr(mask==0)=nan;

smbfinnocorr=0.0;
for i = 252:263
smbfinnocorr = smbfinnocorr + ncread(['hirham.smb.jan.1990.jan.2011.nocorrec000',num2str(i),'.2d.nc'],'smb');
end
smbfinnocorr = rot90(smbfinnocorr);
smbfinnocorr(mask==0)=nan;

% figure,imagesc(smbininocorr),colormap(jet),colorbar;
% figure,imagesc(smbfinnocorr),colormap(jet),colorbar;
data_smb(:,:,1) = smbini;
data_smb(:,:,2) = smbfin;
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
smb_name = {'elevcorr','noelevcorr'};
beta_name = {'celmer_linear'};
smb_legend_name = {'1995elevcorr','2011elevcorr'};
temp_name = {'temp95'};
E_beta_ini_name = {'nothermal_ebet'};
%E_path_name = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/E_results';
E_path_name = 'J:\ASF_data\E_results';

%*****************************configure***********************************%
  %grey=[0.8,0.8,0.8]; 
  grey=[0.6,0.6,0.6]; 
 
  %linecolor= ['c', 'm','b','g'];
  x_B3 = [322 352 377 413 419 426 430];
  y_B3 = [188 248 265 270 268 265 259];
  x_lab_B3 = ((340:20:420)-322)*0.4;
  xlimit_B3 = [322 430];
  datalimit_B3 = [-130 800]; 
  %datalimit_B3 = [0 70]; 
  %datalimit_B3 = [-200 800]; 
  betalimit_B3 = [-0.4 0.8];
  xrange = 'x_B3'; yrange = 'y_B3';
  xlimit = 'xlimit_B3'; datalimit = 'datalimit_B3';
  x_lab = 'x_lab_B3';betalimit = 'betalimit_B3';
  %legend_posi =  [0.51,0.6,0.1,0.3];
  color_data = ['y','r'];
  legend_posi = 'SouthWestOutside';
  counter = 1;
  [ArrayXG(1).data,ArrayYG(1).data,ArrayCG(1).data]=improfile(topg,eval(xrange),...
     eval(yrange),'bicubic');
  [ArrayXG(2).data,ArrayYG(2).data,ArrayCG(2).data]=improfile(topg+thk,eval(xrange),...
     eval(yrange),'bicubic');
  figure,%subplot(2,1,i)
  fillX = [ArrayXG(1).data; flipud(ArrayXG(1).data)];
  fillY = [ArrayCG(2).data; flipud(ArrayCG(1).data)];
  Hfill(1) = fill(fillX,fillY,grey);hold on;
if 0  
%*****************************BISICLES************************************%
for i = 1: length(beta_name)
for j = 1: length(smb_name)
%*******************************Elmer**************************************%

filename = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,cell2mat(E_beta_ini_name(i)),'.ep.mat'];
filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'_ini.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_','20thsmb.ep10.mat'];
load(cell2mat(filename));
Eini=load(cell2mat(filename_ini));

 
Elmer_xvel_bed = matrixdata(2).data_layer(1).data;
Elmer_yvel_bed = matrixdata(3).data_layer(1).data;
Elmer_zvel_bed = matrixdata(4).data_layer(1).data;
Elmer_xvel_surf = matrixdata(2).data_layer(2).data;
Elmer_yvel_surf = matrixdata(3).data_layer(2).data;
Elmer_zvel_surf = matrixdata(4).data_layer(2).data;
Elmer_FS_surf = matrixdata(9).data_layer(1).data;
% data_Ebeta = matrixdata(8).data_layer(1).data;
% Elmer_beta = data_Ebeta;
% data_Ebeta = log10(data_Ebeta);
Elmer_xvel_ini_bed = Eini.matrixdata(2).data_layer(1).data;
Elmer_yvel_ini_bed = Eini.matrixdata(3).data_layer(1).data;
Elmer_zvel_ini_bed = Eini.matrixdata(4).data_layer(1).data;
Elmer_xvel_ini_surf = Eini.matrixdata(2).data_layer(2).data;
Elmer_yvel_ini_surf = Eini.matrixdata(3).data_layer(2).data;
Elmer_zvel_ini_surf = Eini.matrixdata(4).data_layer(2).data;
Elmer_FS_ini_surf = Eini.matrixdata(9).data_layer(1).data;

Elmer_FS_ini_surf (mask==0) = NaN;
Elmer_FS_surf (mask==0) = NaN;
Elmer_FS_change = Elmer_FS_surf - Elmer_FS_ini_surf;

bVel_Elmer = sqrt(Elmer_xvel_bed.^2 +Elmer_yvel_bed.^2+Elmer_zvel_bed.^2);
bVel_Elmer (mask==0) = NaN;

bVel_ini_Elmer = sqrt(Elmer_xvel_ini_bed.^2 + Elmer_yvel_ini_bed.^2+Elmer_zvel_ini_bed.^2);
bVel_ini_Elmer (mask==0) = NaN;
Elmer_bVel_change = bVel_Elmer - bVel_ini_Elmer;

fVel_Elmer = sqrt(Elmer_xvel_surf.^2 +Elmer_yvel_surf.^2+Elmer_zvel_surf.^2);
fVel_Elmer (mask==0) = NaN;

fVel_ini_Elmer = sqrt(Elmer_xvel_ini_surf.^2 + Elmer_yvel_ini_surf.^2+Elmer_zvel_ini_surf.^2);
fVel_ini_Elmer (mask==0) = NaN;
Elmer_fVel_change = fVel_Elmer - fVel_ini_Elmer;
    %*************************************************************************
  fprintf('plotting profile...\n');
  
 [ArrayX(1).data,ArrayY(1).data,ArrayC(1).data]=improfile(Elmer_FS_surf,eval(xrange),...
     eval(yrange),'bicubic');
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data]=improfile(Elmer_fVel_change,eval(xrange),...
     eval(yrange),'bicubic');
 [ArrayX(3).data,ArrayY(3).data,ArrayC(3).data]=improfile(data_smb(:,:,j),eval(xrange),...
     eval(yrange),'bicubic');
 
  Hfill(counter + 1) = plot(ArrayX(1).data,ArrayC(1).data);hold on;
  set(Hfill(counter + 1),'LineWidth',2,'Color',color_data(j));
  xlabel('Distance along section (km)','FontSize',15);
  
     
    [AX,Hplot1(counter).data,Hplot2(counter).data] = plotyy(ArrayX(2).data,ArrayC(2).data,ArrayX(3).data,ArrayC(3).data,'plot');
    hold on;
    smb_legend_name_input = cell2mat(smb_legend_name(j));
    
    legend_FS(counter) = {[' FS_{',cell2mat(smb_name(j)),'}']};
    legend_vel(counter) = {[' \Deltau_{',cell2mat(smb_name(j)),'}']};
    legend_smb(counter) = {[' SMB_{',smb_legend_name_input,'}']};
    
    
    set(AX(1),'YLim',eval(datalimit),'XLim',eval(xlimit),'YColor','k',...
     'XTickLabel',eval(x_lab),'FontSize',15,'ytick',linspace(datalimit_B3(1),datalimit_B3(2),7));
    set(get(AX(1),'Ylabel'),'string','Elevtion (m), or Speed (m/a)','FontSize',15); 
    %ylabel('Elevtion (m), or Speed (m/a)');
    
    set(Hplot1(counter).data,'Color',get(Hfill(counter+1),'Color'),'LineStyle', '-','LineWidth',1.5);
    
    
    set(AX(2),'YLim',eval(betalimit),'XLim',eval(xlimit),'YColor','k',...
      'XTickLabel',eval(x_lab),'FontSize',15,'ytick',linspace(betalimit_B3(1),betalimit_B3(2),5));
    
    set(get(AX(2),'Ylabel'),'string','SMB (m.i.eq)','FontSize',15); 
    set(Hplot2(counter).data,'Color',get(Hfill(counter+1),'Color'),'LineStyle', '--','LineWidth',1.5);
    
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','b','LineStyle','-.','LineWidth',1.5);hold on;
     %ylabel('\beta (MPa/a/m)');
   %set(gca,'XTickLabel',eval(x_lab),'FontSize',15);
   counter = counter + 1;
end
end
 if 0
    LH=legend([Hfill(1);Hfill(2);Hfill(3);...
        Hplot1(1).data;Hplot1(2).data;...
        Hplot2(1).data;Hplot2(2).data;...
        lineH],...
        'Glacier_{initial}',...
        cell2mat(legend_FS(1)),cell2mat(legend_FS(2)),...
        cell2mat(legend_vel(1)),cell2mat(legend_vel(2)),...
        cell2mat(legend_smb(1)),cell2mat(legend_smb(2)),...
        'water line',...
        'Location',legend_posi);
    set(LH,'color','w','FontSize',15);
 end
    
%      'Glacier_{initial}','FS_E_l_m_e_r_/_I_c_e','FS_B_I_S_I_C_L_E_S',...
%      '\Deltau_E_l_m_e_r_/_I_c_e','lg(\beta_E_l_m_e_r_/_I_c_e)','\Deltau_B_I_S_I_C_L_E_S',' ','water line',...
%      'Location',legend_posi);
%    set(LH,'color','none')
% data_B = BISICLES_beta(mask==2);
% data_E = Elmer_beta(mask==2);
% % data_B = fVel_BISICLES;
% % data_E = fVel_Elmer;
% data_B = data_B.*(10^-6);
% data_B(isnan(data_B))=0;
% data_E(isnan(data_E))=0;
% ErrC =  (sum(abs(data_B-data_E))/sum(abs(data_E)))*100

save('profile_B3.mat','ArrayX','ArrayY');

end
smb_20th = ncread('ASF_SMB90_monthly_extrop.nc','smb_extrop');

nodes =69003;% sum(sum(mask~=0));%69003 %elmer nodes
smblimit =  [-600 300];
xlimit = [0 204];
for i = 60:264
    tmp = ncread(strcat('hirham.smb.jan.1990.jan.2011.elevcorrec',sprintf('%06d',i),'.2d.nc'),'smb');
    data=sum(tmp,1);
    data=sum(data,2);
    sum_smb(i-59)=data;
    tmp = ncread(strcat('hirham.smb.jan.1990.jan.2011.nocorrec',sprintf('%06d',i),'.2d.nc'),'smb');
    data=sum(tmp,1);
    data=sum(data,2);
    sum_nocorr_smb(i-59)=data;
    smb_HH = load(['hirham_smb_jan_1990_jan_2011_',num2str(i),'.txt']);
    sum_smb_HH(i-59)=sum(smb_HH(:,4));
end
smb_20thsum(1,1:length(sum_nocorr_smb))=((sum(sum(smb_20th))*(1000/910.0))./nodes).*1e3;
sum_smb = ((sum_smb.*(1000/910))./nodes).*1e3;
sum_nocorr_smb = ((sum_nocorr_smb.*(1000/910.0))./nodes).*1e3;
sum_smb_HH = sum_smb_HH/length(smb_HH(:,4));
figure,
plot(0:length(sum_smb)-1,sum_smb,'.-r','MarkerSize',20,'LineWidth',1.5);hold on;
%plot(0:length(sum_smb)-1,sum_smb,'-k','LineWidth',1.5);hold on;
plot(0:length(sum_nocorr_smb)-1,sum_nocorr_smb,'.-b','MarkerSize',20,'LineWidth',1.5);hold on;
%plot(0:length(sum_nocorr_smb)-1,sum_nocorr_smb,'-m','LineWidth',1.5);hold on;
plot(0:length(sum_smb_HH)-1,sum_smb_HH,'-sk','MarkerSize',15,'LineWidth',1.5);hold on;
plot(0:length(sum_nocorr_smb)-1,smb_20thsum,'--','LineWidth',2,'Color',grey);hold on;
set(gca,'YLim',smblimit,'XLim',xlimit,'xtick',[0,48,96,144,192],...
             'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],...
     'FontSize',40,'ytick',[-600,-300,0,300]);
 ylabel('Mean Nodal SMB(mm w.e. a^{-1})');
 grid on;
legend('SMB_{elevcorr}','SMB_{noelevcorr}','SMB_{HH}','SMB_{90s}','Location','SouthEast','Location','SouthEast');
if 0
r_elecorr = abs((sum_smb-sum_smb_HH)./sum_smb_HH);
figure,plot(1:length(r_elecorr),r_elecorr);
r_elecorr_max = max(r_elecorr(r_elecorr<1));
figure,plot(1:length(r_elecorr(r_elecorr<1)),r_elecorr(r_elecorr<1));

r_noelecorr = abs((sum_nocorr_smb-sum_smb_HH)./sum_smb_HH);
figure,plot(1:length(r_noelecorr),r_noelecorr);
r_noelecorr_max = max(r_noelecorr(r_noelecorr<1));
figure,plot(1:length(r_noelecorr(r_noelecorr<1)),r_noelecorr(r_noelecorr<1));
end

