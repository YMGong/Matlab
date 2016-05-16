   
%*************************************************************************%
%****************** Comparision between Elmer/Ice and BISICLES ************
%*************************************************************************%
clf
 close all;
% axi_units
 clear all;

%*************************************************************************%
fprintf('data loading...\n');
grey = [0.5,0.5,0.5];
path3=('/wrk/ygong/ASFforwardFiles/Invers/ASF_Adjoint_Elmer_400/');
% path2=('E:\ASF\Results_comp\results\Elmer_cost_800\');
% path3=('E:\ASF\Results_comp\results\BISICELS_cost_200\');
path4=('/wrk/ygong/BISICLES/ASF_simulation/ASF_BISICLES_400/');

path2 = 'J:\ASF_data\B_results\';
path1 = '../';

save_path = 'J:\ASF_data\model_compr\';
load profile_B3.mat;
load ASF_B3_mask.mat;  
load ([path1,'ASF_Elmer_400_t6_t5.mat']);
load ([path1,'asf_adjoint_ti_2011_6.ep.mat']);
%load ([path3,'BISICLES_200_temp_t4_C.mat']);
%load ([path4,'BISICLES_400_temp_t4_C.mat']);

N_coordi = load('profileCoordi_y_N_.mat');
B3_coordi = load('profileCoordi_y_B3_.mat');

%BISICLES_95Beta_t4_C.beta = ncread([path2,'BISICLES_400_temp_velomatch_t4_C.nc'],'basal_friction');
%BISICLES_2011Beta_t4_C.beta = ncread([path2,'BISICLES_400_temp_velomatch_t4_vel2011_C.nc'],'basal_friction');
BISICLES_95Beta_t4_C.beta = ncread([path2,'BISICLES_400_Beta_t4_new.nc'],'C');
BISICLES_2011Beta_t4_C.beta = ncread([path2,'BISICLES_400_Beta_t4_vel2011_new.nc'],'C');

%BISICLES_95Beta_t4_C.xvel = ncread([path2,'BISICLES_400_temp_velomatch_t4_C.nc'],'xfVel');
%BISICLES_2011Beta_t4_C.xvel = ncread([path2,'BISICLES_400_temp_velomatch_t4_vel2011_C.nc'],'xfVel');
BISICLES_95Beta_t4_C.xvel = ncread([path2,'BISICLES_400_Beta_t4_new.nc'],'xVels');
BISICLES_2011Beta_t4_C.xvel = ncread([path2,'BISICLES_400_Beta_t4_vel2011_new.nc'],'xVels');

%BISICLES_95Beta_t4_C.yvel = ncread([path2,'BISICLES_400_temp_velomatch_t4_C.nc'],'yfVel');
%BISICLES_2011Beta_t4_C.yvel = ncread([path2,'BISICLES_400_temp_velomatch_t4_vel2011_C.nc'],'yfVel');
BISICLES_95Beta_t4_C.yvel = ncread([path2,'BISICLES_400_Beta_t4_new.nc'],'yVels');
BISICLES_2011Beta_t4_C.yvel = ncread([path2,'BISICLES_400_Beta_t4_vel2011_new.nc'],'yVels');

%BISICLES_95Beta_temp_t4_C.temp = ncread([path2,'BISICLES_400_temp_velomatch_t4_C.nc'],'temperatureBase';
%BISICLES_2011Beta_temp_t4_C.temp = ncread([path2,'BISICLES_400_temp_velomatch_t4_vel2011_C.nc'],'temperatureBase');
BISICLES_95Beta_temp_t4_C.temp = ncread([path2,'BISICLES_400_temp_t3_new.nc'],'temperatureBase');
BISICLES_2011Beta_temp_t4_C.temp = ncread([path2,'BISICLES_400_temp_t3_vel2011_new.nc'],'temperatureBase');
%---------------------------------
Elmer_2011Beta_t6_beta = matrixdata(1).data_layer(1).data;
Elmer_2011Beta_t6_velu_surface = matrixdata(2).data_layer(11).data;
Elmer_2011Beta_t6_velv_surface = matrixdata(3).data_layer(11).data;
Elmer_2011Beta_t6_velw_surface = matrixdata(4).data_layer(11).data;

%---------------------------------
Elmer_2011Beta_t5_temp_bed = matrixdata(5).data_layer(11).data;
%---------------------------------

obvel95x = ncread([path2,'ASF95Efor_BISICLES.nc'],'obvelx');
obvel95y = ncread([path2,'ASF95Efor_BISICLES.nc'],'obvely');
obvel2011x = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'xvel');
obvel2011x = rot90(obvel2011x);
obvel2011y = ncread([path2,'ASF2011Rotatfor_BISICLES.nc'],'yvel');
obvel2011y = rot90(obvel2011y);

velo95 = sqrt(obvel95x.^2 + obvel95y.^2);
velo95 = rot90(velo95);
velo95(mask == 0)= NaN;
%velo95(velo95<3) = NaN;

velo2011 = sqrt(obvel2011x.^2 + obvel2011y.^2);
velo2011(mask == 0)= NaN;
%velo2011(velo2011<3) = NaN;

obvel = {'velo95','velo2011'};

costB_basename = 'cost_velomatch_';
costE_basename = 'Cost_asf_adjoint_ti_';

tempB_basename = {'BISICLES_95Beta_temp_t4_C.temp','BISICLES_2011Beta_temp_t4_C.temp'};
tempE_basename = {'Elmer_400_t5_temp_bed','Elmer_2011Beta_t5_temp_bed'};

betaB_basename = {'BISICLES_95Beta_t4_C.beta','BISICLES_2011Beta_t4_C.beta'};
betaE_basename = {'Elmer_400_t6_beta','Elmer_2011Beta_t6_beta'};

velB_basename = {'BISICLES_95Beta_t4_C','BISICLES_2011Beta_t4_C'};
velB_compon = {'.xvel','.yvel'};

velE_basename = {'Elmer_400_t6_vel','Elmer_2011Beta_t6_vel'};
velE_compon = {'u_surface','v_surface','w_surface'};
% set utm coordi
xs = 545;  xe = 755;
ys = 8400; ye = 8950;
stepD = 0.4;row_exp=384; colume_exp=640;

Xq  = linspace(xs,xe,floor((xe-xs)/stepD+1));
Yq  = linspace(ys,ye,floor((ye-ys)/stepD+1));
Yq(end:-1:1) = Yq(1:end);

y_utm=Yq;y_utm(end:row_exp) = Yq(end,end):-stepD:(Yq(end,end)-(row_exp-size(Yq,2))*stepD);
x_utm=rot90(Xq,-1);x_utm(end:colume_exp) = Xq(end,end):stepD:(Xq(end,end)+(colume_exp-size(Xq,2))*stepD);
fprintf('cost function... \n');
%x_lab = 'total number of iterations';
%y_lab = 'J^o';
file_num = 4;
if 0
for i = 1:file_num
    Cost_BISICLES_95beta(i).data = load([path4,costB_basename,'t',num2str(i),'_vel95_new_C.txt']);
    Cost_BISICLES_2011beta(i).data = load([path4,costB_basename,'t',num2str(i),'_new_C.txt']);
    Cost_Elmer_95beta(i).data = load([path3,costE_basename,num2str(2*(i-1)),'.sif.dat']);
    Cost_Elmer_2011beta(i).data = load([path3,costE_basename,'2011_',num2str(2*(i-1)),'.sif.dat']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Preparing data... \n');
% plot the misfit not the totel cost function
%pdfure,
iter_pre = 0;iter=0;
for i=1:4
    iter = size(Cost_BISICLES_2011beta(i).data(:,1))+iter;
    h1=plot((1+iter_pre):iter,Cost_BISICLES_2011beta(i).data(:,3),'s','MarkerSize',5,'LineWidth',1.5);%,'Color',[abs((4-i)/3) abs((3-i)/3) abs((2-i)/3)]);
    hold on;
    h2=plot((1+iter_pre):iter,Cost_BISICLES_95beta(i).data(:,3),'o','MarkerSize',5,'LineWidth',1.5);%,'Color',[abs((4-i)/3) abs((3-i)/3) abs((2-i)/3)]);
    hold on;
    %LH = legend('BISICLES\_\beta_2_0_1_1','BISICLES\_\beta_9_5');
    
    cost_B_95Beta((1+iter_pre):iter,1) = (1+iter_pre):iter;
    cost_B_95Beta((1+iter_pre):iter,2) = Cost_BISICLES_95beta(i).data(:,3);
    cost_B_2011Beta((1+iter_pre):iter,1) = (1+iter_pre):iter;
    cost_B_2011Beta((1+iter_pre):iter,2) = Cost_BISICLES_2011beta(i).data(:,3);
    iter_pre = size(Cost_BISICLES_2011beta(i).data(:,1))+iter_pre;
end
 grid on;xlabel('Total number of iterations','FontSize',15);ylabel('J^o','FontSize',15);
 %title('cost function BISICLES 400','FontSize',15);
 %set(LH, 'FontSize',15);
 set(gca,'FontSize',15);
 %ylim([3e11 5e12]);
 hold on;
 %saveas(gcf,'.\pics\Cost-function-BISICLES-beta95-beta2011.tif','tif');
 %clf;

%figure;
iter_pre = 0;iter = 0;
for i=1:4
    iter = size(Cost_Elmer_2011beta(i).data(:,1),1)+iter;
    h3=plot((1+iter_pre):iter,Cost_Elmer_2011beta(i).data(:,3),'s','MarkerSize',5,'LineWidth',1.5,'MarkerFaceColor','k');%,'Color',[abs((4-i)/3) abs((3-i)/3) abs((2-i)/3)]);
    hold on;
    h4=plot((1+iter_pre):iter,Cost_Elmer_95beta(i).data(:,3),'o','MarkerSize',5,'LineWidth',1.5,'MarkerFaceColor','k');%,'Color',[abs((4-i)/3) abs((3-i)/3) abs((2-i)/3)]);
    hold on;
    %LH = legend('Elmer/Ice\_\beta_2_0_1_1','Elmer/Ice\_\beta_9_5');
    cost_E_95Beta((1+iter_pre):iter,1) = (1+iter_pre):iter;
    cost_E_95Beta((1+iter_pre):iter,2) = Cost_Elmer_95beta(i).data(:,3);
    cost_E_2011Beta((1+iter_pre):iter,1) = (1+iter_pre):iter;
    cost_E_2011Beta((1+iter_pre):iter,2) = Cost_Elmer_2011beta(i).data(:,3);
    iter_pre = size(Cost_Elmer_2011beta(i).data(:,1),1)+iter_pre;
end
 grid on;xlabel('Total number of iterations','FontSize',15);ylabel('J^o','FontSize',15);
 
 %LH = legend('BISICLES\_\beta_2_0_1_1','BISICLES\_\beta_9_5',...
   %  'Elmer/Ice\_\beta_2_0_1_1','Elmer/Ice\_\beta_9_5');
 %set(LH, 'FontSize',15);
 h=[h1,h2,h3,h4];
 legend(h,'BISICLES\_\beta_2_0_1_1','BISICLES\_\beta_9_5','Elmer/Ice\_\beta_2_0_1_1','Elmer/Ice\_\beta_9_5');
 set(gca,'FontSize',15);
 %title('cost function Elmer ');
 ylim([3e11 5e12]);

 saveas(gcf,'.\pics\Cost_function_ElmerIce_BISICLE_beta95.tif','tif');
end

 clf;
 
 for  timeperiod = 1%2:-1:1
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
boxcolor = {'k','y','r','m'}; gray = [0.6 0.6 0.6];
%*************************************************************************%


fprintf('velocity ...\n');
c_units = 'm a^{-1}';colorbar_lab = (0:1)*100;

switch timeperiod 
    case 1
    colorbarlim = [0 20];
    otherwise
    colorbarlim = [0 50];
end

gridon = 1; lineon = 1;boxon = 1;
%BISICLES data needs to be rotated
if 0
    velo95(velo95<3) = NaN;
    velo2011(velo95<3) = NaN;
for ii = 240:270
    for jj = 175:200
        if velo95(ii,jj)<16
            velo95(ii,jj)=NaN;
        end
    end
end
velo95(300:end,280:500)=NaN;
for ii = 270:300
    for jj = 370:500
        if velo95(ii,jj)<25
            velo95(ii,jj)=NaN;
        end
    end
end

for ii = 430:490
    for jj = 80:200
        if velo95(jj,ii)<16
            velo95(jj,ii)=NaN;
        end
    end
end
end

%2 for 2011beta, 1 for 95beta
%velB_basename = {'BISICLES_200_Beta_t4_C','BISICLES_400_Beta_t4_C'};
%velB_compon = {'.xvel','.yvel'};
velB_basename = {'BISICLES_95Beta_t4_C','BISICLES_2011Beta_t4_C'};
velB_compon = {'.xvel','.yvel'};
data_B = ((eval(cell2mat([velB_basename(timeperiod),velB_compon(1)]))).^2 + ...
         (eval(cell2mat([velB_basename(timeperiod),velB_compon(2)]))).^2).^(1/2); %- velo;

%data_B = (abs(data_B - velo)./abs(velo))*100; %compute relative the difference
%#####################
data_B = rot90(data_B);
data_B = abs(data_B - (eval(cell2mat(obvel(timeperiod)))));
data_B(mask==0)=NaN;

data_E = ((eval([cell2mat(velE_basename(timeperiod)),cell2mat(velE_compon(1))])).^2 +...
         (eval([cell2mat(velE_basename(timeperiod)),cell2mat(velE_compon(2))])).^2 +...
         (eval([cell2mat(velE_basename(timeperiod)),cell2mat(velE_compon(3))])).^2).^(1/2);
%data_E = (abs(data_E - velo)./abs(velo))*100; %compute the relative difference
data_E = abs(data_E - (eval(cell2mat(obvel(timeperiod)))));
data_E(mask==0)=NaN;

data_diff = data_B-data_E;

if 0
for i = 1:length(xlabs)
figure,H=imagesc(data_E);colorbar;hold on;
set(H,'alphadata',~isnan(data_E));
if strcmp(cell2mat(xlabs(i)),'x_lab')
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
    
end 
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
edittingASFplot(eval(cell2mat(xaxis(i))),eval(cell2mat(yaxis(i))),colorbarlim,...
    eval(cell2mat(xlabs(i))),eval(cell2mat(ylabs(i))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(i)));
saveas(gcf,['.\pics\',cell2mat(velE_basename(timeperiod)),'_velabsdiff_',num2str(i),'.tif'],'tif');

end


for j = 1:length(xlabs)
figure,H=imagesc(data_B);colorbar;hold on;
set(H,'alphadata',~isnan(data_E));
if strcmp(cell2mat(xlabs(j)),'x_lab')
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
 
end
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
edittingASFplot(eval(cell2mat(xaxis(j))),eval(cell2mat(yaxis(j))),colorbarlim,...
    eval(cell2mat(xlabs(j))),eval(cell2mat(ylabs(j))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(j)));
saveas(gcf,['.\pics\',cell2mat(velB_basename(timeperiod)),'_velabsdiff_',num2str(j),'.tif'],'tif');
i=i+1;
end

end
%************************************************************************
fprintf('velocity diff...\n');
%*************************************************************************%
xaxi = [150 500];yaxi = [50 370];
%x_lab = (200:100:500)*0.400; y_lab = (50:50:370)*0.400;
x_tick = linspace(xaxi(1),xaxi(2),6);
y_tick = linspace(yaxi(1),yaxi(2),5);
x_lab = round(x_utm(x_tick));
y_lab = round(y_utm(y_tick));
gridon = 1; lineon = 0;boxon = 1;

switch timeperiod 
    case 1
    colorbarlim_diff = [-4 4];
    otherwise
    colorbarlim_diff = [-15 15];
end
for k = 1:length(xlabs)
colorbarlim=[min(min(eval(cell2mat(obvel(timeperiod))))),max(max(eval(cell2mat(obvel(timeperiod)))))-60];
if strcmp(cell2mat(xlabs(k)),'x_lab')
    figure,H=imagesc(eval(cell2mat(obvel(timeperiod))));
    colorbar;colormap(jet);hold on;set(H,'alphadata',~isnan(eval(cell2mat(obvel(timeperiod)))));
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','Color',grey,'LineWidth',2),hold on;
    plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);hold on;
%     plot(b3boxx,linspace(b3boxy(1),b3boxy(1),11),'--k','LineWidth',2);hold on;
%     plot(b3boxx,linspace(b3boxy(end),b3boxy(end),11),'--k','LineWidth',2);hold on;
%     plot(linspace(b3boxx(1),b3boxx(1),11),b3boxy,'--k','LineWidth',2);hold on;
%     plot(linspace(b3boxx(end),b3boxx(end),11),b3boxy,'--k','LineWidth',2);hold on;
    %plot(N_coordi.ArrayX(1).data,N_coordi.ArrayY(1).data,'LineStyle','--','Color',gray,'LineWidth',2);
    %plot(B3_coordi.ArrayX(1).data,B3_coordi.ArrayY(1).data,'LineStyle',':','Color',gray,'LineWidth',2);
    edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),colorbarlim,...
    eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
    x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(k)));
    hx=xlabel('Easting UTM33X(km)','FontSize',20);%,'Position',[0.5,-0.01]);
    posx=get(hx,'pos');
    set(hx,'pos',posx+[0.0 -14.5 0]);
    hy=ylabel('Northing UTM33X(km)','FontSize',20);%,'Position',[-0.17,0.1]);
    posy=get(hy,'pos');
    set(hy,'pos',posy+[-38.0 0 0]); 
else
    if 0
figure,H=imagesc(data_diff);colorbar;hold on;set(H,'alphadata',~isnan(data_diff));
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),colorbarlim_diff,...
    eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(k)));

colormap(darkb2r(-10,10));
    end
end

%saveas(gcf,['.\pics\',cell2mat(betaE_basename(timeperiod)),'_veldiffB-E_',num2str(k),'.tif'],'tif');
i=i+1;

end
if 0
figure,H=imagesc(data_diff);colorbar;hold on;set(H,'alphadata',~isnan(data_B));
plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
plot(N_coordi.ArrayX(1).data,N_coordi.ArrayY(1).data,'LineStyle','--','Color',gray,'LineWidth',2);
plot(B3_coordi.ArrayX(1).data,B3_coordi.ArrayY(1).data,'LineStyle',':','Color',gray,'LineWidth',2);
edittingASFplot(eval(cell2mat(xaxis(1))),eval(cell2mat(yaxis(1))),colorbarlim_diff,...
    eval(cell2mat(xlabs(1))),eval(cell2mat(ylabs(1))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(1)));
saveas(gcf,['.\pics\',cell2mat(betaB_basename(timeperiod)),'_veldiffB-E_',num2str(1),'.tif'],'tif');

% data_B = ((eval(cell2mat([velB_basename(timeperiod),velB_compon(1)]))).^2 + ...
%           (eval(cell2mat([velB_basename(timeperiod),velB_compon(2)]))).^2).^(1/2); %- velo;
% data_B = rot90(data_B);
% data_B(mask==0)=NaN;
% 
% data_E = ((eval([cell2mat(velE_basename(timeperiod)),cell2mat(velE_compon(1))])).^2 +...
%          (eval([cell2mat(velE_basename(timeperiod)),cell2mat(velE_compon(2))])).^2 +...
%          (eval([cell2mat(velE_basename(timeperiod)),cell2mat(velE_compon(3))])).^2).^(1/2);
% data_E(mask==0)=NaN;
%  figure(i),imagesc(data_B),colorbar;title('BISICLES');
%  figure(i+1),imagesc(data_E),colorbar;title('ELmer/Ice');
%  figure(i+2),imagesc(velo),colorbar;title('observation')

%******************************************************************************************8
end
% zoom in!!!
xaxi = [150 500];yaxi = [50 370];
x_tick = linspace(xaxi(1),xaxi(2),6);
y_tick = linspace(yaxi(1),yaxi(2),5);
x_lab = round(x_utm(x_tick));
y_lab = round(y_utm(y_tick));
%x_lab = (200:100:500)*0.400; y_lab = (50:50:370)*0.400;
b3boxx = linspace(305,445,11);b3boxy = linspace(180, 320,11);
fprintf('beta ...\n');
clear data_B;
clear data_E;
clear data_diff;
%*************************************************************************%
c_units = 'MPa a^{-1} m^{-1}';colorbar_lab = -4:0;
gridon = 1; lineon = 0;boxon = 1;
switch timeperiod 
    case 1
    colorbarlim = [-6 0];
    otherwise
    colorbarlim = [-6 0];
end
%x_lab = ('80', '120', '160', ' '); y_lab = (50:50:384)*0.4; 
%BISICLES data needs to be rotated

%2 for 400m, 1 for the other  log10(Bate_2011.*(10^-6))
data_B = log10(eval(cell2mat(betaB_basename(timeperiod))).*(10^-6));
data_B = rot90(data_B);
data_B(mask==0)=NaN;

data_E = eval(cell2mat(betaE_basename(timeperiod)));
data_E(mask==0)=NaN;

data_diff = data_B - data_E;
data_diff(mask == 0)=NaN;


for k = 1:length(xlabs)

if strcmp(cell2mat(xlabs(k)),'x_lab')
    figure,H=imagesc(data_E);colorbar;colormap(jet);hold on;set(H,'alphadata',~isnan(data_E));
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','Color',grey,'LineWidth',2),hold on;
    plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);hold on;
%     plot(b3boxx,linspace(b3boxy(1),b3boxy(1),11),'--k','LineWidth',2);hold on;
%     plot(b3boxx,linspace(b3boxy(end),b3boxy(end),11),'--k','LineWidth',2);hold on;
%     plot(linspace(b3boxx(1),b3boxx(1),11),b3boxy,'--k','LineWidth',2);hold on;
%     plot(linspace(b3boxx(end),b3boxx(end),11),b3boxy,'--k','LineWidth',2);hold on;
    %plot(N_coordi.ArrayX(1).data,N_coordi.ArrayY(1).data,'LineStyle','--','Color',gray,'LineWidth',2);
    %plot(B3_coordi.ArrayX(1).data,B3_coordi.ArrayY(1).data,'LineStyle',':','Color',gray,'LineWidth',2);
    edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),colorbarlim,...
    eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
    x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(k)));
    hx=xlabel('Easting UTM33X(km)','FontSize',20);%,'Position',[0.5,-0.01]);
    posx=get(hx,'pos');
    set(hx,'pos',posx+[0.0 -14.5 0]);
    hy=ylabel('Northing UTM33X(km)','FontSize',20);%,'Position',[-0.17,0.1]);
    posy=get(hy,'pos');
    set(hy,'pos',posy+[-37.0 0 0]); 
    
else
if 0
figure,H=imagesc(data_diff);colorbar;hold on;set(H,'alphadata',~isnan(data_diff));
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
edittingASFplot(eval(cell2mat(xaxis(k))),eval(cell2mat(yaxis(k))),[-1 1],...
    eval(cell2mat(xlabs(k))),eval(cell2mat(ylabs(k))),...
    x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(k)));

colormap(darkb2r(-1,1));
end
end
end

saveas(gcf,['.\pics\',cell2mat(betaE_basename(timeperiod)),'_lg(beta)_',num2str(k),'.tif'],'tif');


if 0
figure,H=imagesc(data_B);colorbar;hold on;set(H,'alphadata',~isnan(data_B));
plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','Color',grey,'LineWidth',2),hold on;
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);
%plot(N_coordi.ArrayX(1).data,N_coordi.ArrayY(1).data,'LineStyle','--','Color',gray,'LineWidth',2);
%plot(B3_coordi.ArrayX(1).data,B3_coordi.ArrayY(1).data,'LineStyle',':','Color',gray,'LineWidth',2);
edittingASFplot(eval(cell2mat(xaxis(1))),eval(cell2mat(yaxis(1))),colorbarlim,...
    eval(cell2mat(xlabs(1))),eval(cell2mat(ylabs(1))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(1)));
xlabel('km','FontSize',20,'Position',[0.5,-0.06]);
saveas(gcf,['.\pics\',cell2mat(betaB_basename(timeperiod)),'_lg(beta)_',num2str(1),'.tif'],'tif');
end
fprintf('temperature ...\n');
clear data_B;
clear data_E;
clear data_diff;
%*************************************************************************%
 %tempB_basename = {'BISICLES_200_temp_t4_C(12).temp','BISICLES_400_temp_t4_C(12).temp'};
 %tempE_basename = {'Elmer_800_t5_temp_bed','Elmer_400_t5_temp_bed'};

c_units = '^oC';colorbar_lab = -10:0;colorbarlim = [-10 0];
gridon = 1; lineon = 0;boxon = 1;

data_B = eval(cell2mat(tempB_basename(timeperiod)))- 273.15;
data_B = rot90(data_B);
data_B(mask==0)=NaN;

data_E = eval(cell2mat(tempE_basename(timeperiod)))-273.15;
data_E(mask==0)=NaN;

data_diff = data_B - data_E;
data_diff(mask == 0)=NaN;
for p = 1:length(xlabs)

if strcmp(cell2mat(xlabs(p)),'x_lab')
    figure,H=imagesc(data_E);colorbar;colormap(jet);hold on;set(H,'alphadata',~isnan(data_E));
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','Color',grey,'LineWidth',2),hold on;
    plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);hold on;
%     plot(b3boxx,linspace(b3boxy(1),b3boxy(1),11),'--k','LineWidth',2);hold on;
%     plot(b3boxx,linspace(b3boxy(end),b3boxy(end),11),'--k','LineWidth',2);hold on;
%     plot(linspace(b3boxx(1),b3boxx(1),11),b3boxy,'--k','LineWidth',2);hold on;
%     plot(linspace(b3boxx(end),b3boxx(end),11),b3boxy,'--k','LineWidth',2);hold on;
    %plot(N_coordi.ArrayX(1).data,N_coordi.ArrayY(1).data,'LineStyle','--','Color',gray,'LineWidth',2);
    %plot(B3_coordi.ArrayX(1).data,B3_coordi.ArrayY(1).data,'LineStyle',':','Color',gray,'LineWidth',2);
    edittingASFplot(eval(cell2mat(xaxis(p))),eval(cell2mat(yaxis(p))),colorbarlim,...
    eval(cell2mat(xlabs(p))),eval(cell2mat(ylabs(p))),...
    x_tick,y_tick,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(p)));
    hx=xlabel('Easting UTM33X(km)','FontSize',20);%,'Position',[0.5,-0.01]);
    posx=get(hx,'pos');
    set(hx,'pos',posx+[0.0 -14.5 0]);
    hy=ylabel('Northing UTM33X(km)','FontSize',20);%,'Position',[-0.17,0.1]);
    posy=get(hy,'pos');
    set(hy,'pos',posy+[-37.0 0 0]); 
    
else
    if 0
    figure,H=imagesc(data_diff);colorbar;hold on;set(H,'alphadata',~isnan(data_diff));
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    edittingASFplot(eval(cell2mat(xaxis(p))),eval(cell2mat(yaxis(p))),[-2 2],...
    eval(cell2mat(xlabs(p))),eval(cell2mat(ylabs(p))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(p)));

colormap(darkb2r(-2,2));
    end
end

saveas(gcf,['.\pics\',cell2mat(tempE_basename(timeperiod)),'_temphomo_',num2str(p),'.tif'],'tif');
i=i+1;
end
if 0
figure,H=imagesc(data_B);colorbar;hold on;set(H,'alphadata',~isnan(data_B));

plot(B3_FoPr(1,:),B3_FoPr(2,:),'-g','LineWidth',2),hold on;
plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
edittingASFplot(eval(cell2mat(xaxis(1))),eval(cell2mat(yaxis(1))),colorbarlim,...
    eval(cell2mat(xlabs(1))),eval(cell2mat(ylabs(1))),...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,cell2mat(boxcolor(1)));
saveas(gcf,['.\pics\',cell2mat(tempB_basename(timeperiod)),'_temphomo_',num2str(1),'.tif'],'tif');

 end
 end

%load ASF_B3_mask.mat;
load '../asf_b3_mask.txt';

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
colorbarlim = [-6,0];
colorbar_lab = (0:1)*100;
grey = [0.6 0.6 0.6];
%c_units = 'm/yr';
c_units = 'MPa a^{-1} m^{-1}';
Tinterv = 16.0;
betaE_basename = {'Elmer_400_t6_beta','Elmer_2011Beta_t6_beta'};
slope =  (10.^Elmer_2011Beta_t6_beta - 10.^Elmer_400_t6_beta) / Tinterv;
Fric = max(10.^-30,10.^Elmer_2011Beta_t6_beta + slope * 1.0);
load ASF_B3_mask.mat;
Fric(mask~=2)=10.^Elmer_2011Beta_t6_beta(mask~=2);
betaextrop = log10(Fric);
fprintf('Plotting... \n');

figure;
H=imagesc(betaextrop);colorbar;%title('Elmer: Vel_{surf}-Vel_{surf-control}');
colormap(jet);
hold on;

set(H,'alphadata',mask~=0);
    plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
    plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
    plot(I_xllarge:435,linspace(I_yularge,I_yularge,size(I_xllarge:435,2)),'-.m','LineWidth',2),hold on;
    plot(linspace(I_xllarge,I_xllarge,size(I_yularge:309,2)),I_yularge:309,'-.m','LineWidth',2),hold on;
    plot(I_xl:427,linspace(I_yu,I_yu,size(I_xl:427,2)),'-.b','LineWidth',2),hold on;
    plot(linspace(I_xl,I_xl,size(I_yu:310,2)),I_yu:310,'-.b','LineWidth',2),hold on;
    %plot(ArrayX(1).data,ArrayY(1).data,'--w','LineWidth',2);
    edittingASFplot(xaxi,yaxi,colorbarlim(1,:),xlab,ylab,0,0,...
    colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');

%figure,imagesc(betaextrop),colormap(jet),colorbar;


 
 
 
 