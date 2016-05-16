clear all;
close all;
load smb0304_seasonal.txt;
load('aust_250.txt');
season = {'winter','summer'};
c_units = 'mm a^{-1}';
xaxi = [605 745];yaxi = [8802 8930];
x_tick = linspace(xaxi(1),xaxi(2),6);
y_tick = linspace(yaxi(1),yaxi(2),5);
Cstep = 6;
grey = [0.6 0.6 0.6];
samppoint = [705610,8839203];
x_utm = smb0304_seasonal(:,1)/1000;
y_utm = smb0304_seasonal(:,2)/1000;
xutm_Th = x_utm(~isnan(smb0304_seasonal(:,3)));
yutm_Th = y_utm(~isnan(smb0304_seasonal(:,3)));
DIM=2;

aust_250(aust_250==-9999)=nan;
elev_Th = reshape(flipud(aust_250),length(smb0304_seasonal),1);
elev_Th = elev_Th(~isnan(smb0304_seasonal(:,3)));
%figure,imagesc(aust_250);colormap(jet),colorbar;
% creat a colormap for slope
rgb = [ ...
    %94   79   162
    0     0     0
   153   153   153
   50    136   189
   0     255   255
   230   245   152
   254   224   139
   253   174    97
   244   109    67
   213    62    79
   158     1    66  ] / 255;
for i= 1:size(season,2)
%load elev_HH_interpTh.mat;
corr=load(['HH_Th_correlation_',cell2mat(season(i)),'_rCorr.mat']);
SMB = smb0304_seasonal(:,DIM+i);
Clabel=round(linspace(min(SMB),max(SMB),Cstep));
figure,
opengl software
PlotScatterPoints(x_utm,y_utm,SMB,'jet',Cstep,c_units,Clabel,6);%set(gca,'color',grey);
cbh = findobj( 0, 'tag', 'Colorbar' );
hold on;

SMB2d = flipud(reshape(SMB,518,538));

SMB = SMB(~isnan(SMB)); 

Isamppoint = find(corr.yutm_HH==samppoint(2));
Redius = corr.Redius_HH(Isamppoint)/1000;

coord_Th(:,1) = xutm_Th - corr.crd_HH(Isamppoint).coor(1)/1000;
coord_Th(:,2) = yutm_Th - corr.crd_HH(Isamppoint).coor(2)/1000;
dist = (coord_Th(:,1).^2 + coord_Th(:,2).^2).^(1/2);
pairs_smb = SMB((dist<Redius)|(dist==Redius));
pairs_elev = elev_Th((dist<Redius)|(dist==Redius));
plot(xutm_Th((dist<Redius)|(dist==Redius)),yutm_Th((dist<Redius)|(dist==Redius)),...
    '.','MarkerFaceColor',grey,'MarkerEdgeColor',grey,'markersize',2);hold on;

set(gca,'Xlim',xaxi,'Ylim',yaxi,...
    'xtick',x_tick,'ytick',y_tick,...
     'FontSize',15);
 xlabh = xlabel('Easting UTM33X (km)');
 ylabh = ylabel('Northing UTM33X (km)');
 %xlabh = xlabel(axi_units);
 set(xlabh,'FontSize',15);
 set(ylabh,'FontSize',15);
 grid on;
%  plot good and bad points in b&w
% for j=1:size(corr.slope_HH,1)
%     if corr.slope_HH(j) ==0.0
%         plot(corr.xutm_HH(j)/1000,corr.yutm_HH(j)/1000,'.k','markersize',10),hold on;
%     else
%         plot(corr.xutm_HH(j)/1000,corr.yutm_HH(j)/1000,'.w','markersize',10),hold on;
%     end
% end
% 
% PlotScatterPoints(corr.xutm_HH/1000,corr.yutm_HH/1000,corr.slope_HH,...
%     darkb2r(min(corr.slope_HH),max(corr.slope_HH)),...
%     0,'','',10);

% cheated for winter slope in order to set bed points to be 0;
%there are two points that have negtive value

PlotScatterPoints(corr.xutm_HH/1000,corr.yutm_HH/1000,max(0,corr.slope_HH),...
    rgb,0,'','',10);hold on;
% the colormap has changed for some reason, now change it back
negtiveI = find(corr.slope_HH<0);
plot(corr.xutm_HH(negtiveI)/1000,corr.yutm_HH(negtiveI)/1000,'w.','MarkerSize',10);

%set(gca,'color',grey);

% this is just for plotting puposes; I did
% 'regression(pairs_elev,pairs_smb,'one');'for the really downscaling!!!
[r,slope,offset]=regression(pairs_smb,pairs_elev,'one');

figure,plot(pairs_smb,pairs_elev,'.','MarkerFaceColor',grey,'MarkerEdgeColor',grey);hold on;
plot(pairs_smb,slope*pairs_smb+offset,'-r');
text(min(pairs_smb),min(pairs_elev)+5,...
    ['$$\mathsf{ {1 \over k} =',num2str(slope),';r^2=',num2str(r^2), '}$$'],...
    'interpreter','latex','FontSize',15);
set(gca,'FontSize',15);
 xlabh = xlabel('SMB (mm)');
 ylabh = ylabel('Elevation (m)');
 %xlabh = xlabel(axi_units);
 set(xlabh,'FontSize',15);
 set(ylabh,'FontSize',15);
% corr.slope_HH(Isamppoint)
end



figure,
PlotScatterPoints(corr.xutm_HH/1000,corr.yutm_HH/1000,max(0,corr.slope_HH),...
    rgb,...
    8,'(m/mm)',linspace(0,max(corr.slope_HH),8),10);
% PlotScatterPoints(corr.xutm_HH/1000,corr.yutm_HH/1000,corr.slope_HH,...
%     rgb,...%darkb2r(min(corr.slope_HH),max(corr.slope_HH)),...
%     8,'(m/mm)',linspace(min(corr.slope_HH),max(corr.slope_HH),8),10);
