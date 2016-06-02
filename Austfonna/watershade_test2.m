%************************************************************************%
%1.You need to have the GIS tool box, which can be find in :
%/home/ygong/Documents/Acamedia/Coding/matlab_exe;
%make sure that you add the dir properly
%2.You need to creat a file which can be read by GRIDobj like the way below,
%which has been commented out;
%3.You need to emerge small drainage basins together by hand to creat the one you want;
%4.The color code of the first figure is not the index of the grid point, but you can randomly change one color, 
%plot the field again and use the index in that figure.
%02/06/2016, Yongmei Gong
%************************************************************************%

clear all;
close all;

 addpath D:\'Program Files'\Dropbox\matlab\topotoolbox-2.0-r
 addpath D:\'Program Files'\Dropbox\matlab\topotoolbox-2.0-r\topoapp\
 addpath D:\'Program Files'\Dropbox\matlab\topotoolbox-2.0-r\tools_and_more\
 fprintf('loading data...\n');
 load('MaskASF_Elmer.mat');mask = Mask_ASF_Elmer;
 x_utm=ncread('ASF95Efor_BISICLES.nc','x');
 y_utm=ncread('ASF95Efor_BISICLES.nc','y');
%fprintf('creating the dem file for GRIDobj.../n')
% I=rot90(topo+thk);
% %I=BISICLES_freesurface_surf;
% I(I<0)=0;
% I(mask==0)=-9999;
% 
% %I(isnan(I))=0;
% % 
% % Info=imfinfo(BISICLES_freesurface_surf);
% % 
% % if Info.BitDepth>8   
% fid=fopen('ASF_DEM.txt','w');
% fprintf(fid,['ncols    640\n','nrows   384\n','xllcorner   545000.0000\n','yllcorner   8800000.0000\n'....
%         'cellsize   400\n','NODARA_value   -9999\n']);
% for i=1:size(I,1)
%     fprintf(fid,[num2str(I(i,:)),'\n'],'%5.5f');
% end
% fclose(fid);
dem = GRIDobj('ASF_DEM.txt');
hillshade(dem)
FD  = FLOWobj(dem,'preprocess','c');
D   = drainagebasins(FD);
%     imageschs(dem,D);
%******
fprintf('Searching for B3...\n');
B3=D.Z;

     B3(D.Z==8)=300;%
     B3(D.Z==17)=300;%
     B3(D.Z==16)=300;%
     B3(D.Z==13)=300;%
     B3(D.Z==12)=300;%
     B3(D.Z==138)=300;%
     B3(D.Z==119)=300;%
     B3(D.Z==95)=300;%
     B3(D.Z==48)=300;%
     B3(D.Z==72)=300;%
     B3(D.Z==55)=300;%
    B3(D.Z==44)=300;%
    B3(D.Z==33)=300;%
    B3(D.Z==29)=300;%
    B3(D.Z==27)=300;%
    B3(D.Z==25)=300;%
    B3(D.Z==23)=300;%
    B3(D.Z==21)=300;%
    B3(D.Z==22)=300;%
 
figure,imagesc(B3),colorbar;

mask(B3==300)=2;
figure(3),imagesc(mask);hold on;

contour_index=contourc(mask,2); 
ASF_index = contour_index(:,2:contour_index(2,1)+1);
B3_index = contour_index(:,contour_index(2,1)+3:end);
 plot(ASF_index(1,:),ASF_index(2,:),'k','LineWidth',2);
% plot(B3_index(1,:),B3_index(2,:),'b','LineWidth',2);

B3_FoPr=round(B3_index);

index_x=B3_FoPr(1,:)>396;
index_y=B3_FoPr(2,:)>252;

index_B3_marine=index_x&index_y;

B3_marine_FoPr=B3_FoPr(:,index_B3_marine);

B3_side_FoPr=B3_FoPr(:,~index_B3_marine);

B3_side_utm_e(:,1)=x_utm(B3_side_FoPr(1,:),:);
B3_side_utm_e(:,2)=y_utm(B3_side_FoPr(2,:),:);
size_BS=size(B3_side_utm_e,1);
B3_side_utm = B3_side_utm_e(1:size_BS-1,:);
size(B3_side_utm)
[B3_side_utm,IC1,IA1] = unique(B3_side_utm,'rows','stable');
size(B3_side_utm)

B3_marine_utm(:,1)=x_utm(B3_marine_FoPr(1,:),:);
B3_marine_utm(:,2)=y_utm(B3_marine_FoPr(2,:),:);
size(B3_marine_utm)
[B3_marine_utm,IC,IA]=unique(B3_marine_utm,'rows','stable');
size(B3_marine_utm)

plot(B3_side_FoPr(1,:),B3_side_FoPr(2,:),'b.');hold on;
plot(B3_marine_FoPr(1,:),B3_marine_FoPr(2,:),'r.');hold on;
figure(4),plot(B3_side_utm(:,1),B3_side_utm(:,2),'.r','LineWidth',2);
%*******
fprintf('saving...\n');
filename = 'ASF_B3_mask.mat';
save(filename,'mask','ASF_index','B3_side_FoPr','B3_marine_FoPr','B3_FoPr');
save('B3_footprint.mat','B3_side_utm','B3_marine_utm');




