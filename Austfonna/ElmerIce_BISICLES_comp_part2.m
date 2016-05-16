%*************************************************************************%
%****************** Comparision between Elmer/Ice and BISICLES ************
%****************** basal velocity profile and everything *****************
%*************************************************************************%

clf;
close all;
clear all;

%*************************************************************************%
fprintf('data loading...\n');
addpath 'E:\ASF\Results_comp\results\Elmer_linearbeta\';
save_path=('E:\ASF\Results_comp\results\pics\');
load ASF_B3_mask.mat;
topg = ncread('ASF95Efor_BISICLES.nc','topg');
topg = rot90(topg);
thk = ncread('ASF95Efor_BISICLES.nc','thk');
thk = rot90(thk);
load ('asf_transient_linearbeta_abs_1990_2011.ep.mat');
beta95=load ('ASF_Elmer_400_t6_t5.mat');
beta2011=load('asf_transient_2011beta_abs_1990_2011.ep.mat');

fprintf('velocity ...\n');

%2 for 400m, 1 for the other
data_B = (beta95.Elmer_400_t6_velu_bed.^2+...
    beta95.Elmer_400_t6_velv_bed.^2+...
    beta95.Elmer_400_t6_velw_bed.^2).^(1/2); %- velo 95;
data_B(mask==0)=NaN;

data_E = (beta2011.matrixdata(2).data_layer(1).data.^2 +...
         beta2011.matrixdata(3).data_layer(1).data.^2 +...
         beta2011.matrixdata(4).data_layer(1).data.^2).^(1/2);
data_E(mask==0)=NaN;

data_Bbeta = beta95.Elmer_400_t6_beta;
data_Bbeta(mask==0)=NaN;

data_Ebeta = beta2011.matrixdata(9).data_layer(1).data;
data_Ebeta(mask==0)=NaN;

%*************************************************************************
  fprintf('plotting profile...\n');
  grey=[0.8,0.8,0.8]; linecolor= ['c', 'm'];x_lab_N = (1:5:31)*0.4;x_lab_B3 = (3:10:43)*0.4;
  %x_B3 = [390 398 405 412 421 427 431]; y_B3 = [267 268 269 268 264 258 253];
  %x_N = [377 403 424 435 452]; y_N = [162 137 104 89 69];
  x_N = [414 449]; y_N = [167 67];
  x_B3 = [387 436];y_B3 = [275 261];
  xlimit_N = [414 446];
  xlimit_B3 = [387 432];
  ylimit_N = [-120 600]; betalimit_N = [-5 0];
  ylimit_B3 = [-150 650]; betalimit_B3 = [-9.5 -2];
  xrange = {'x_N', 'x_B3'}; yrange = {'y_N', 'y_B3'};
  xlimit = {'xlimit_N', 'xlimit_B3'}; ylimit = {'ylimit_N', 'ylimit_B3'};
  betalimit = {'betalimit_N', 'betalimit_B3'};
  x_lab = {'x_lab_N','x_lab_B3'};
  for i = 2 %B3 1for N
 [ArrayX(1).data,ArrayY(1).data,ArrayC(1).data]=improfile(topg,eval(cell2mat(xrange(i))),...
     eval(cell2mat(yrange(i))),'bicubic');
 [ArrayX(2).data,ArrayY(2).data,ArrayC(2).data]=improfile(topg+thk,eval(cell2mat(xrange(i))),...
     eval(cell2mat(yrange(i))),'bicubic');
 [ArrayX(3).data,ArrayY(3).data,ArrayC(3).data]=improfile(data_B,eval(cell2mat(xrange(i))),...
     eval(cell2mat(yrange(i))),'bicubic');
 [ArrayX(4).data,ArrayY(4).data,ArrayC(4).data]=improfile(data_E,eval(cell2mat(xrange(i))),...
     eval(cell2mat(yrange(i))),'bicubic');
 [ArrayX(5).data,ArrayY(5).data,ArrayC(5).data]=improfile(data_Bbeta,eval(cell2mat(xrange(i))),...
     eval(cell2mat(yrange(i))),'bicubic');
 [ArrayX(6).data,ArrayY(6).data,ArrayC(6).data]=improfile(data_Ebeta,eval(cell2mat(xrange(i))),...
     eval(cell2mat(yrange(i))),'bicubic');
  figure(i),%subplot(2,1,i)
  fillX = [ArrayX(1).data; flipud(ArrayX(1).data)];
  fillY = [ArrayC(2).data; flipud(ArrayC(1).data)];
  Hfill = fill(fillX,fillY,grey);hold on;
  set(Hfill,'LineWidth',2);
  xlabel('Distance along section (km)','FontSize',15);
  legend(Hfill,'glacier');
 %set(gca,'XTickLabel',eval(cell2mat(x_lab(i))),'FontSize',15);
  for j = 1:2
  [AX,Hplot1(j).data,Hplot2(j).data] = plotyy(ArrayX(j+2).data,ArrayC(j+2).data,ArrayX(j+4).data,ArrayC(j+4).data,'plot');hold on;
   set(AX(1),'YLim',eval(cell2mat(ylimit(i))),'XLim',eval(cell2mat(xlimit(i))),'YColor','k',...
       'XTickLabel',eval(cell2mat(x_lab(i))),'FontSize',15);
   ylabel('Elevtion (m), or Speed (m/a)');
   set(Hplot1(j).data,'Color',linecolor(j),'LineStyle', '-','LineWidth',1.5);
   set(AX(2),'YLim',eval(cell2mat(betalimit(i))),'XLim',eval(cell2mat(xlimit(i))),'YColor','k',...
       'XTickLabel',eval(cell2mat(x_lab(i))),'FontSize',15);
   if j==2
   else
   set(get(AX(2),'Ylabel'),'string','lg(\beta) (MPa/a/m)','FontSize',15) 
   end
   set(Hplot2(j).data,'Color',linecolor(j),'LineStyle', '--','LineWidth',1.5);
   end
   lineH=line('XData',eval(cell2mat(xlimit(i))),'YData',[0 0],'Color','b','LineStyle','-.');
   set(gca,'XTickLabel',eval(cell2mat(x_lab(i))),'FontSize',15);
   legend([Hfill;Hplot1(1).data;Hplot2(1).data;Hplot1(2).data;Hplot2(2).data;lineH],...
      'glacier','u_b_9_5','lg(\beta_9_5)','u_b_2_0_1_1','lg(\beta_2_0_1_1)','water line','Location','NorthWest');
   %legend(lineH,'water line');
%    saveas(gcf,['.\pics\',cell2mat(xrange(i)),'_profile.fig'],'fig');
%    save(['profileCoordi_',cell2mat(yrange(i)),'_.mat'],'ArrayX','ArrayY');
  end
 if 0
%**************************************************************************
figure(10),imagesc(data_E);
load ASF_B3_mask.mat;
base_name_Elmer = 'asf_transient_nothermal_notrimmed_';
model_name = {'B','E'};
color_name = ['r','b','k','g'];

%*************************************************************************%
xaxi = [1 640]; yaxi = [1 384]; 
xaxi_N = [300 475];yaxi_N = [55 140];
xaxi_S = [170 250];yaxi_S = [225 310];
xaxi_W = [345 470];yaxi_W = [175 280];

x_lab = (100:100:640)*0.4; y_lab = (50:50:384)*0.4; 
x_lab_N = (300:50:475)*0.4; y_lab_N = (60:10:140)*0.4; 
x_lab_S = (180:20:250)*0.4; y_lab_S = (230:10:310)*0.4; 
x_lab_W = (350:20:470)*0.4; y_lab_W = (180:20:280)*0.4; 
%let the lab of the axis to have the right scale 
axi_units = ' ';
xaxis = {'xaxi','xaxi_N','xaxi_S','xaxi_W'};
yaxis = {'yaxi','yaxi_N','yaxi_S','yaxi_W'};
xlabs = {'x_lab','x_lab_N','x_lab_S','x_lab_W'};
ylabs = {'y_lab','y_lab_N','y_lab_S','y_lab_W'};
boxcolor = {'k','y','r','m'}; gray = [0.8 0.8 0.8];
%*************************************************************************%

fprintf('velocity ...\n');
c_units = '(m/a)';colorbar_lab = 0:200;colorbarlim = [0 200];
gridon = 1; lineon = 1;boxon = 1;
flag = 1;color_flag = 1;
pic = 1;

path_name = ['Elmer_',cell2mat(model_name(2)),'_Beta_',cell2mat(model_name(2)),'_temp'];
load (['.\',path_name,'\',base_name_Elmer,...
                lower(cell2mat(model_name(2))),'b',lower(cell2mat(model_name(2))),'t',...
                '_1990_2011.ep.mat'],'matrixdata');
            xvel = matrixdata(2).data_layer(11).data;
            yvel = matrixdata(3).data_layer(11).data;
            zvel = matrixdata(4).data_layer(11).data;
            vel = ((xvel).^2 +(yvel).^2 +(zvel).^2).^(1/2);
            vel(mask==0)=NaN;

for i = 1:length(model_name)
   for j = 1:length(model_name)
      B_path_name = ['BISICLES_',cell2mat(model_name(i)),'_Beta_',cell2mat(model_name(j)),'_temp']; 
      E_path_name = ['Elmer_',cell2mat(model_name(i)),'_Beta_',cell2mat(model_name(j)),'_temp']; 
      if strcmp(B_path_name,'BISICLES_E_Beta_B_temp')
           % do nothing       
      else
        B_xvel = ncread(['.\',B_path_name,'\',cell2mat(model_name(i)),'_Beta_',cell2mat(model_name(j)),'_Temp','.nc'],'xfVel');
        B_yvel = ncread(['.\',B_path_name,'\',cell2mat(model_name(i)),'_Beta_',cell2mat(model_name(j)),'_Temp','.nc'],'yfVel');
        data_B = ((B_xvel).^2 + (B_xvel).^2).^(1/2); 
        data_B = rot90(data_B);
        data_B(mask==0)=NaN;
        
        data_diff = data_B - vel;
        for p = 1:length(xlabs)
         
        if strcmp(cell2mat(xlabs(p)),'x_lab')
            figure(pic),imagesc(data_B),colorbar;hold on;
            plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',3),hold on;
            plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',3),hold on;
            edittingASFplot(eval(cell2mat(xaxis(p))),eval(cell2mat(yaxis(p))),colorbarlim,...
            eval(cell2mat(xlabs(p))),eval(cell2mat(ylabs(p))),...
            colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(p)));
        else
            figure(pic),imagesc(data_diff),colorbar;hold on;
            plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',3),hold on;
            edittingASFplot(eval(cell2mat(xaxis(p))),eval(cell2mat(yaxis(p))),[-30 30],...
            eval(cell2mat(xlabs(p))),eval(cell2mat(ylabs(p))),...
            colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(p)));
            colormap(darkb2r(-30,30));
        end
            saveas(gcf,['.\pics\',B_path_name,'_vel99to11_',num2str(pic),'.fig'],'fig');
            pic = pic +1;% loop over plots
        end
        
       end
      if strcmp(E_path_name,'Elmer_B_Beta_E_temp')
          
           load (['.\',E_path_name,'\',base_name_Elmer,...
                lower(cell2mat(model_name(i))),'b',lower(cell2mat(model_name(j))),'t',...
                '_1990_2011.ep.mat'],'matrixdata');
            E_xvel = matrixdata(2).data_layer(11).data;
            E_yvel = matrixdata(3).data_layer(11).data;
            E_zvel = matrixdata(4).data_layer(11).data;
      
        data_E = ((E_xvel).^2 +(E_yvel).^2 +(E_zvel).^2).^(1/2);
        data_E(mask==0)=NaN;
        data_diff = data_E - vel;
        
        for k = 1:length(xlabs)
        if strcmp(cell2mat(xlabs(k)),'x_lab')
            figure(pic),imagesc(data_E),colorbar;hold on;    
            plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',3),hold on;
            plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',3),hold on;
            edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),colorbarlim,...
                eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
            colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(p)));
        else
            figure(pic),imagesc(data_diff),colorbar;hold on; 
            plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',3),hold on;
            edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),[-30 30],...
                eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
            colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(k)));
            colormap(darkb2r(-30,30));
            
        end
        saveas(gcf,['.\pics\',E_path_name,'_vel99to11_',num2str(pic),'.fig'],'fig');
        pic = pic +1;% loop over plots
        end
         
      elseif strcmp(E_path_name,'Elmer_E_Beta_E_temp')
         load (['.\',E_path_name,'\',base_name_Elmer,...
                lower(cell2mat(model_name(i))),'b',lower(cell2mat(model_name(j))),'t',...
                '_1990_2011.ep.mat'],'matrixdata');
            E_xvel = matrixdata(2).data_layer(11).data;
            E_yvel = matrixdata(3).data_layer(11).data;
            E_zvel = matrixdata(4).data_layer(11).data;
      
        data_E = ((E_xvel).^2 +(E_yvel).^2 +(E_zvel).^2).^(1/2);
        data_E(mask==0)=NaN;
               
        for k = 1:length(xlabs)
         figure(pic),imagesc(data_E),colorbar;hold on;      
        if strcmp(cell2mat(xlabs(k)),'x_lab')
         plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',3),hold on;
    
        end
        plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',3),hold on;
        edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),colorbarlim,...
            eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
        colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(k)));
        saveas(gcf,['.\pics\',E_path_name,'_vel99to11_',num2str(pic),'.fig'],'fig');
        pic = pic +1;% loop over plots
        end 
      end  
        
   end
end

end







