%************************************************************************%
clf
clear all;
close all;
fprintf('data loading...\n');
load ASF_B3_mask.mat;
N_coordi = load('profileCoordi_y_N_.mat');
B3_coordi = load('profileCoordi_y_B3_.mat');
topg = ncread('J:\ASF_data\B_results\ASF95Efor_BISICLES.nc','topg');
thk = ncread('J:\ASF_data\B_results\ASF95Efor_BISICLES.nc','thk');
 surf = topg + thk;
 surf = rot90(surf);
 surf(mask == 0)= NaN;
 topg=rot90(topg);
 topg(mask == 0)=NaN;

fprintf('parameterization...\n');
%N_line = 15;%the number of the coutour line you want to drew
spaceSurf = 50;
spaceTopg = 20;
lineSpaceSurf = (min(min(surf)):spaceSurf:max(max(surf)));
lineSpaceTopg = (min(min(topg)):spaceTopg:max(max(topg)));
%xaxi = [305 445];yaxi = [180 320]; % B3

%graph edit
tpg = surf;
flag = 1;
xaxi = [150 500];yaxi = [50 370]; %ASF
%xaxi = [305 445];yaxi = [180 320];%B3
xs = 545;  xe = 755;
ys = 8400; ye = 8950;
stepD = 0.4;row_exp=384; colume_exp=640;
Xq  = linspace(xs,xe,floor((xe-xs)/stepD+1));
Yq  = linspace(ys,ye,floor((ye-ys)/stepD+1));
Yq(end:-1:1) = Yq(1:end);
y_utm=Yq;y_utm(end:row_exp) = Yq(end,end):-stepD:(Yq(end,end)-(row_exp-size(Yq,2))*stepD);
x_utm=rot90(Xq,-1);x_utm(end:colume_exp) = Xq(end,end):stepD:(Xq(end,end)+(colume_exp-size(Xq,2))*stepD);
colorbarlim=[min(min(tpg)),max(max(tpg))];

x_tick = linspace(xaxi(1),xaxi(2),6);
y_tick = linspace(yaxi(1),yaxi(2),5);
xlab = round(x_utm(x_tick));
ylab = round(y_utm(y_tick));
gridon = 1; lineon = 0;boxon = 1;
%colorbarlim=[0,5000;0,15000;-10,750];
%colorbarlim = [0,5000;0,15000];
%colorbar_lab = (0:1)*100;
grey = [0.5 0.5 0.5];
colorbar_lab = 0:100:800;

c_units = 'm';axi_units = ' ';
b3boxx = linspace(305,445,11);b3boxy = linspace(180, 320,11);

fprintf('plot and add contour line...\n');  

  figure,H = imagesc(tpg);set(H,'alphadata',~isnan(tpg));hold on;
  %hAx1 = axes('Parent',figure,'Color','none');
  plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','Color',grey,'LineWidth',2),hold on;
  plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
  xlim(xaxi);ylim(yaxi);
  
if flag  
 
 [c,h] = contour(tpg,lineSpaceSurf,'k','LineWidth',1.2);axis ij; 
 texth=clabel(c,h,'manual','fontsize', 10); 
 for i=1:size(texth) 
     textstr=get(texth(i),'String');
     textnum=str2num(textstr); 
     textstrnew=sprintf('%0.0f', textnum);
     set(texth(i),'String',textstrnew);
 end
 plot(b3boxx,linspace(b3boxy(1),b3boxy(1),11),'--k','LineWidth',2);hold on;
 plot(b3boxx,linspace(b3boxy(end),b3boxy(end),11),'--k','LineWidth',2);hold on;
 plot(linspace(b3boxx(1),b3boxx(1),11),b3boxy,'--k','LineWidth',2);hold on;
 plot(linspace(b3boxx(end),b3boxx(end),11),b3boxy,'--k','LineWidth',2);hold on;
else 
    colorbarlim = [-150, 300];
    [c1,h1] = contour(surf,lineSpaceSurf,'w','LineWidth',1.2);axis ij; 
    [c,h] = contour(tpg,lineSpaceTopg,'k','LineWidth',1.2);axis ij; 
end
 plot(ASF_index(1,:),ASF_index(2,:),'-k','LineWidth',2),hold on;
 plot(B3_FoPr(1,:),B3_FoPr(2,:),'-','LineWidth',2,'Color',grey),hold on;
 edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
 x_tick,y_tick,...
 colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
 hx=xlabel('Easting UTM33X(km)','FontSize',20);%,'Position',[0.5,-0.01]);
 posx=get(hx,'pos');
 set(hx,'pos',posx+[0.0 -14.5 0]);
 %set(hx,'pos',posx+[0.0 -5 0]);
 hy=ylabel('Northing UTM33X(km)','FontSize',20);%,'Position',[-0.17,0.1]);
 posy=get(hy,'pos');
 set(hy,'pos',posy+[-37.0 0 0]); 
 
 figure,H = imagesc(tpg);set(H,'alphadata',~isnan(tpg));hold on;
 edittingASFplot(xaxi,yaxi,colorbarlim,xlab,ylab,...
 x_tick,y_tick,...
 colorbar_lab,c_units,axi_units,gridon,lineon,boxon,'k');
 %hx=xlabel('Easting UTM33X(km)','FontSize',20);%,'Position',[0.5,-0.01]);
 %posx=get(hx,'pos');
 %set(hy,'pos',posy+[-15.0 0 0]); 
%  Bed rock