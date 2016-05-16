%close all;
clear all;
fprintf('Loading data... \n');
%path2 = '/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/';
path2 = 'J:\ASF_data\B_results\';
%addpath /wrk/ygong/BISICLES/ASF_simulation/PostProcessing/data4bisicles/;
addpath J:\ASF_data\B_results
load ASF_B3_mask.mat;
topg = ncread([path2,'ASF95Efor_BISICLES.nc'],'topg');
topg = rot90(topg);
thk = ncread([path2,'ASF95Efor_BISICLES.nc'],'thk');
thk = rot90(thk);
% smbini_data = smbini.matrixdata(1).data_layer.data/910.0;
% smbini_data(mask==0)=nan;
% smbfin_data = smbfin.matrixdata(1).data_layer.data/910.0;
% smbfin_data(mask==0)=nan;
if 0
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
end
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
smb_legend_name = {'elevcorr\_1995','elevcorr\_2011'};
temp_name = {'temp95'};
E_beta_ini_name = {'nothermal_ebet'};
%E_path_name = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/';
%B_path_name = '/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/data4bisicles/';
B_path_name = 'J:\ASF_data\B_results\';
%*****************************configure***********************************%
  grey=[0.8,0.8,0.8]; 
  Azero = 273.15; 
  %linecolor= ['c', 'm','b','g'];
  x_B3 = [322 352 377 413 419 426 430];
  y_B3 = [188 248 265 270 268 265 259];
  x_lab_B3 = ((340:20:420)-322)*0.4;
  xlimit_B3 = [322 430];
  datalimit_B3 = [-130 800]; 
  colorbarlim = [262 273] - Azero;
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
  layer = 11;
  
  [ArrayXG(1).data,ArrayYG(1).data,ArrayCG(1).data]=improfile(topg,eval(xrange),...
     eval(yrange),'bicubic');
  [ArrayXG(2).data,ArrayYG(2).data,ArrayCG(2).data]=improfile(topg+thk,eval(xrange),...
     eval(yrange),'bicubic');
  figure(1),%subplot(2,1,i)
  fillX = [ArrayXG(1).data; flipud(ArrayXG(1).data)];
  fillY = [ArrayCG(2).data; flipud(ArrayCG(1).data)];
  Hfill(1) = fill(fillX,fillY,grey);hold on;
  set(Hfill(1),'LineWidth',1.5);

%*****************************BISICLES************************************%
for i = 1: layer

%*******************************Elmer**************************************%
temp_layer = ncread([B_path_name,'frc.velomatch.monthly0000016.elmer.2d.nc'],['temperature',sprintf('%06d',i-1)]);
%temp_layer = ncread([B_path_name,'frc.velomatch.monthly0000000.elmer.2d.nc'],['temperature',sprintf('%06d',i-1)]);
temp_layer = rot90(temp_layer) - Azero;
%*************************************************************************
  fprintf('plotting profile...\n');
  
 [ArrayX(1).data,ArrayY(1).data,ArrayC(1).data]=improfile(temp_layer,eval(xrange),...
     eval(yrange),'bicubic');
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data]=improfile(topg+(layer-i)*(thk/(layer-1)),eval(xrange),...
     eval(yrange),'bicubic');
 
    colormap(jet(256)); % or whatever colormap you want
    ax=surface('XData',  [ArrayX(1).data ArrayX(1).data],'YData',[ArrayC(2).data ArrayC(2).data],...
        'ZData',0*[ArrayX(1).data ArrayX(1).data] ,'CData',[ArrayC(1).data ArrayC(1).data],'EdgeColor','flat','LineWidth',3);
    
    axc=colorbar;
    ylabel(axc,'^oC','FontSize',15,'FontWeight','bold');
    %set(axc, 'XTick', [colorbarlim(1), sum(colorbarlim)/2, colorbarlim(2)])
    
    opengl software
    set(gca,'YLim',eval(datalimit),'XLim',eval(xlimit),'Clim',colorbarlim,...
     'YColor','k','XColor','k','box','on',...
     'XTickLabel',eval(x_lab),'FontSize',15,'FontWeight','bold',...
     'ytick',linspace(datalimit_B3(1),datalimit_B3(2),7));
     set(get(gca,'Ylabel'),'string','Elevtion (m)','FontSize',15,'FontWeight','bold'); 
     set(get(gca,'Xlabel'),'string','Distance (km)','FontSize',15,'FontWeight','bold'); 
     %caxis([263,273]);
    
     hold on;
  
    %ylabel('Elevtion (m), or Speed (m/a)');XTickLabel
    if 0
    set(Hplot1(counter).data,'Color',Hfill(counter+1).Color,'LineStyle', '-','LineWidth',1.5);
    
    
    set(AX(2),'YLim',eval(betalimit),'XLim',eval(xlimit),'YColor','k',...
      'XTickLabel',eval(x_lab),'FontSize',15,'ytick',linspace(betalimit_B3(1),betalimit_B3(2),5));
    
    set(get(AX(2),'Ylabel'),'string','SMB (m.i.eq)','FontSize',15); 
    set(Hplot2(counter).data,'Color',Hfill(counter+1).Color,'LineStyle', '--','LineWidth',1.5);
    
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','b','LineStyle','-.','LineWidth',1.5);hold on;
     %ylabel('\beta (MPa/a/m)');
   %set(gca,'XTickLabel',eval(x_lab),'FontSize',15);
   counter = counter + 1;
    end
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','b','LineStyle','-.','LineWidth',1.5);hold on;
end
legend([Hfill(1);lineH],'Glacier','water line');
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

%save('profile_B3.mat','ArrayX','ArrayY');



