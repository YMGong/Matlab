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
time_name = {'12to13'};
smb_name = {'20thsmb.ep6'};
beta_name = {'celmer_triangelnofrc','celmer_triangellargenofrc'};
beta_legend_name = {'nofrc1','nofrc2'};
temp_name = {'temp95'};
E_beta_ini_name = {'nothermal_ebet'};
E_path_name = 'J:\ASF_data\E_results\';
%figure, imagesc(topg),colormap(jet),colorbar;

%*****************************configure***********************************%
  grey=[0.8,0.8,0.8]; 
  %linecolor= ['c', 'm','b','g'];
  x_B3_N = [392 433];y_B3_N = [273 257];
  x_B3_M = [392 428];y_B3_M =[273 286];
  x_B3_S = [392 419];y_B3_S =[272 302];
    
  %x_B3 = [322 352 377 413 419 426 430];
  %y_B3 = [188 248 265 270 268 265 259];
  x_lab_B3 = ((340:20:420)-322)*0.4;
  %xlimit_B3 = [322 430];
  yy1limit_B3 = [-150 300]; 
  %datalimit_B3 = [0 800]; 
  %datalimit_B3 = [-200 800]; 
  %betalimit_B3 = [-9 0];
  savefilename = {'profile_B3_N.mat','profile_B3_M.mat','profile_B3_S.mat'};
  yy2limit_B3 = [0 1200; 0 8000;0 8000];
  xrange = {'x_B3_N','x_B3_M','x_B3_S'}; yrange = {'y_B3_N','y_B3_M','y_B3_S'};
   yy1limit = 'yy1limit_B3';
  x_lab = 'x_lab_B3';yy2limit = 'yy2limit_B3';
 for k = 1:length(xrange)   
  %legend_posi =  [0.51,0.6,0.1,0.3];
  xlimit = cell2mat(xrange(k));
  legend_posi = 'Best';
  counter = 1;
  [ArrayXG(1).data,ArrayYG(1).data,ArrayCG(1).data]=improfile(topg,eval(cell2mat(xrange(k))),...
     eval(cell2mat(yrange(k))),'bicubic');
  [ArrayXG(2).data,ArrayYG(2).data,ArrayCG(2).data]=improfile(topg+thk,eval(cell2mat(xrange(k))),...
     eval(cell2mat(yrange(k))),'bicubic');
  

%*****************************************************************%
for j = 1: length(smb_name)
 figure;
  fillX = [ArrayXG(1).data; flipud(ArrayXG(1).data)];
  fillY = [ArrayCG(2).data; flipud(ArrayCG(1).data)];
  Hfill(1) = fill(fillX,fillY,grey);hold on;
  set(Hfill(1),'LineWidth',2);   
for i = 1: length(beta_name)
%*********************************************************************%

filename = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_',cell2mat(smb_name(j)),'.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,cell2mat(E_beta_ini_name(i)),'.ep.mat'];
filename_ini = [E_path_name,'/',base_name_Elmer,'_95to11_',temp_name,'_celmer_linear_elevcorr_ini.ep10.mat'];
%filename_ini = [E_path_name,'/',base_name_Elmer,'_',time_name,'_',temp_name,'_',cell2mat(beta_name(i)),'_','20thsmb.ep9.mat'];
load(cell2mat(filename));
Eini=load(cell2mat(filename_ini));

 
Elmer_xvel_bed = matrixdata(2).data_layer(1).data;
Elmer_yvel_bed = matrixdata(3).data_layer(1).data;
Elmer_zvel_bed = matrixdata(4).data_layer(1).data;
Elmer_xvel_surf = matrixdata(2).data_layer(2).data;
Elmer_yvel_surf = matrixdata(3).data_layer(2).data;
Elmer_zvel_surf = matrixdata(4).data_layer(2).data;
Elmer_FS_surf = matrixdata(5).data_layer(1).data;

Elmer_xvel_ini_surf = Eini.matrixdata(2).data_layer(2).data;
Elmer_yvel_ini_surf = Eini.matrixdata(3).data_layer(2).data;
Elmer_zvel_ini_surf = Eini.matrixdata(4).data_layer(2).data;
 Elmer_FS_ini_surf = Eini.matrixdata(9).data_layer(1).data;
% % 
 Elmer_FS_ini_surf (mask==0) = NaN;
Elmer_FS_surf (mask==0) = NaN;
Elmer_FS_change = Elmer_FS_surf;%-Elmer_FS_ini_surf;

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
  
 [ArrayX(1).data,ArrayY(1).data,ArrayC(1).data] = improfile(Elmer_fVel_change,eval(cell2mat(xrange(k))),...
    eval(cell2mat(yrange(k))),'bicubic');
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data] = improfile(Elmer_FS_change,eval(cell2mat(xrange(k))),...
     eval(cell2mat(yrange(k))),'bicubic');
 
 ArrayC(1).data(end)=Elmer_fVel_change(eval([cell2mat(yrange(k)),'(end)']),...
     eval([cell2mat(xrange(k)),'(end)']));
  ArrayC(2).data(end)=Elmer_FS_change(eval([cell2mat(yrange(k)),'(end)']),...
     eval([cell2mat(xrange(k)),'(end)']));

    [AX,Hplot1(counter).data,Hplot2(counter).data] =...
        plotyy(ArrayX(2).data(~isnan(ArrayC(2).data)),ArrayC(2).data(~isnan(ArrayC(2).data)),...
        ArrayX(1).data(~isnan(ArrayC(2).data)),ArrayC(1).data(~isnan(ArrayC(2).data)),'plot');
    hold on;
    beta_legend_name_input = cell2mat(beta_legend_name(i));
    
    %legend_beta(counter) = {[' \beta_{',beta_legend_name_input,'}']};
    
    
    set(AX(1),'YLim',eval(yy1limit),'XLim',eval(xlimit),'YColor','k',...
     'xtick',linspace(eval([xlimit,'(1)']),eval([xlimit,'(end)']),6),...
     'XTickLabel',((linspace(eval([xlimit,'(1)']),eval([xlimit,'(end)']),6))-eval([xlimit,'(1)']))*0.4,...
     'FontSize',15,'ytick',linspace(yy1limit_B3(1),yy1limit_B3(2),11));
    set(get(AX(1),'Ylabel'),'string','Elevtion (m)','FontSize',15); 
    set(get(AX(1),'Xlabel'),'string','Distance (km)','FontSize',15); 
        %ylabel('Elevtion (m), or Speed (m/a)');,'XTickLabel',eval(x_lab),...
    set(AX(2),'YLim',yy2limit_B3(k,:),'XLim',eval(xlimit),'YColor','k',...
      'xtick',linspace(eval([xlimit,'(1)']),eval([xlimit,'(end)']),6),...
      'XTickLabel',((linspace(eval([xlimit,'(1)']),eval([xlimit,'(end)']),6))-eval([xlimit,'(1)']))*0.4,...
      'FontSize',15,'ytick',linspace(-1*yy2limit_B3(k,2),yy2limit_B3(k,2),9));
    set(get(AX(2),'Ylabel'),'string','Speed (m/a)','FontSize',15); 
    set(get(AX(2),'Xlabel'),'string','Distance (km)','FontSize',15); 
%     if i==3||i==4
%         
%         set(Hplot2(counter).data,'Color','none','LineStyle', '--','LineWidth',1.5);
%     else   

        set(Hplot1(counter).data,'LineStyle', '--','LineWidth',2);  
        set(Hplot2(counter).data,'Color',Hplot1(counter).data.Color,'LineStyle', '-','LineWidth',2);
        legend_FS(counter) = {['FS_{',beta_legend_name_input,'}']};
        legend_vel(counter) = {['|u|_{',beta_legend_name_input,'}']};

    counter = counter + 1;
 
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','b','LineStyle','-.','LineWidth',1.5);hold on;
     %ylabel('\beta (MPa/a/m)');
     %set(gca,'XTickLabel',eval(x_lab),'FontSize',15);
     %figure,imagesc(data_Ebeta),colorbar,colormap(jet);
end
end
if 0
    LH=legend([Hfill(1);Hfill(2);Hfill(3);Hfill(4);Hfill(5);...
        Hplot1(1).data;Hplot1(2).data;Hplot1(3).data;Hplot1(4).data;...
        Hplot2(1).data;Hplot2(2).data;Hplot2(3).data;Hplot2(4).data;...
        lineH],...
        'Glacier_{initial}',...
        cell2mat(legend_FS(1)),cell2mat(legend_FS(2)),cell2mat(legend_FS(3)),cell2mat(legend_FS(4)),...
        cell2mat(legend_vel(1)),cell2mat(legend_vel(2)),cell2mat(legend_vel(3)),cell2mat(legend_vel(4)),...
        '\beta_{BISICLES}','\beta_E_l_m_e_r_/_I_c_e','','',...
        'sea level',...
        'Location',legend_posi);
    set(LH,'color','w','FontSize',15);
end
   
   LH=legend(gca(figure(1)),[Hfill;...
       Hplot1(1).data;Hplot1(2).data;...
       Hplot2(1).data;Hplot2(2).data;...
       lineH],...
       'Glacier_{initial}',...
       cell2mat(legend_FS(1)),cell2mat(legend_FS(2)),...
       cell2mat(legend_vel(1)),cell2mat(legend_vel(2)),...
       'sea level',...
       'Location',legend_posi);
    set(LH,'color','none','FontSize',15);
   
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

save(cell2mat(savefilename(k)),'ArrayX','ArrayY');
end

