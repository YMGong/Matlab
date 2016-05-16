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
time_name = {'95to11'};
smb_name = {'elevcorr_ini','elevcorr_yr05','elevcorr_yr10','elevcorr'};
beta_name = {'celmer_linear','celmer_stepc95','celmer_stepc2011'};
beta_legend_name = {'linear','step1995','step2011'};
time_legend_name = {'','2005','2010','2011'};
temp_name = {'temp95'};
E_beta_ini_name = {'nothermal_ebet'};
E_path_name = 'J:\ASF_data\E_results';

%*****************************configure***********************************%
  grey=[0.8,0.8,0.8]; 
  %linecolor= ['c', 'm','b','g'];
  x_B3 = [322 352 377 413 419 426 430];
  y_B3 = [188 248 265 270 268 265 259];
  x_lab_B3 = ((340:20:420)-322)*0.4;
  xlimit_B3 = [322 430];
  yy1limit_B3 = [-130 800]; 
  %datalimit_B3 = [0 800]; 
  %datalimit_B3 = [-200 800]; 
  %betalimit_B3 = [-9 0];
  yy2limit_B3 = [-50 40];
  xrange = 'x_B3'; yrange = 'y_B3';
  xlimit = 'xlimit_B3'; yy1limit = 'yy1limit_B3';
  x_lab = 'x_lab_B3';yy2limit = 'yy2limit_B3';
  color = {'r',[0.6 0.4 0.4],'b'};
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
if i==1
    LC = 'k';
else
    LC= cell2mat(color(j));
end
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
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data]=improfile(Elmer_fVel_change,eval(xrange),...
    eval(yrange),'bicubic');
  
    figure(1);
    H(counter) = plot(ArrayX(2).data,ArrayC(2).data,...
         cell2mat(linestyle(i)),'LineWidth',2,'Color',LC);
     
    hold on;
    beta_legend_name_input = cell2mat(beta_legend_name(j));
    time_legend_name_input = cell2mat(time_legend_name(i));
  
    legend_vel(counter) = {[time_legend_name_input,': |u|_{',beta_legend_name_input,'}']};

    counter = counter + 1;
 
    lineH=line('XData',eval(xlimit),'YData',[0 0],'Color','c','LineStyle','-.','LineWidth',1.5);hold on;
  
end

  if i==1
  else
    figure(2),
    H2(counter) = plot(ArrayX(1).data,ArrayC(1).data,...
         cell2mat(linestyle(i)),'LineWidth',2,'Color',LC);hold on;
  end
end


set(gca,'XTickLabel',eval(x_lab),...
    'XLim',eval(xlimit), 'YLim',eval(yy1limit),...
    'FontSize',15);%somehow the y direction of grid has been reversed
    ylabel('Elevtion (m) or Speed (m a^{-1})');
    ylim(eval(yy1limit));
% axis tight;
    % 'YTick',linspace(yy1limit_B3(1),yy1limit_B3(2),11), ...

   LH=legend(gca(figure(1)),[Hfill;H(1);H(2);H(3);H(4);H(6);H(7);H(8);H(10);H(11);H(12);lineH], ...
       'Glacier_{initial}',...
       '|u|_{initial}',...
       cell2mat(legend_vel(2)),cell2mat(legend_vel(3)),cell2mat(legend_vel(4)),...
       cell2mat(legend_vel(6)),cell2mat(legend_vel(7)),cell2mat(legend_vel(8)),...
       cell2mat(legend_vel(10)),cell2mat(legend_vel(11)),cell2mat(legend_vel(12)),...
       'sea level',...
       'Location',legend_posi);
    set(LH,'color','w','FontSize',30);

set(gca(figure(2)),...
    'XTickLabel',eval(x_lab),...
    'XLim',eval(xlimit),...
    'FontSize',15);%somehow the y direction of grid has been reversed
    ylabel('Elevtion Change (m)');xlabel('Distance(km)','FontSize',15);
 ylim(eval(yy2limit));   
    axis tight;

save('profile_B3.mat','ArrayX','ArrayY');


