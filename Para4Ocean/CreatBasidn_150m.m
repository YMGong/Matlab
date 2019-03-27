%************************************************************************%
%1.You need to have the GIS tool box, which can be find in :
%/Users/yongmei/Documents/Academy/Coding/Matlab
%make sure that you add the dir properly
%1.You need to emerge small drainage basins together by hand to creat the one you want;
%2.The color code of the first figure is not the index of the grid point, but you can randomly change one color, 
%plot the field again and use the index in that figure.
%11/01/2018, Yongmei Gong
%************************************************************************%

clear all;
close all;
addpath(genpath('/Users/yongmei/Documents/Academy/Project/PARAMOUR/data_BedMachine/'));
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_BedMachine/Morlighem17_Datasets/Original_dataset/
addpath /Users/yongmei/Library/'Application Support'/MathWorks/'MATLAB Add-Ons'/Toolboxes/'Arctic Mapping Tools'/'Arctic Mapping Tools'

%%
%load the constants for greenland
config_greenland 
fprintf('Data loading ...\n');
%{
I = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','surface'); I = rot90(I); % the surface elevation
xpsn = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','lon'); xpsn = rot90(xpsn); % the x coordinates
ypsn = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','lat'); ypsn = rot90(ypsn); % the y coordinates
mask = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','mask'); mask = rot90(mask); 
% 0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land
I(mask~=2) = 0.0; 
mask_ice = mask; mask_ice(mask~=2) = 0.0; 
%%
fprintf('Parameterization ...\n');
saveGridfile = 0; % 0 = do not save the GridObj file created by mat2gridobj
resl = 150; %150m grid size
%%
fprintf('Find the basins ...\n');
[x,y] = polarstereo_fwd(ypsn,xpsn,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
% convert lat/lon data for a polar stereographic system to decimal degrees
% on maps with negative numbers (-) for S and W. 
xllcorner = x(end,1); % ideally x coordinate of the Center Point of the lower left Corner Grid Cell
yllcorner = y(end,1); % ideally y coordinate of the Center Point of the lower left Corner Grid Cell
%[xutm,yutm] = ll2utm(x,y,'wgs84',19);
dem = mat2gridobj(I,mask_ice, xllcorner, yllcorner, resl, saveGridfile); %creat a GRIDobj
hillshade(dem)
FD  = FLOWobj(dem,'preprocess','c');
DB   = drainagebasins(FD);
%DB   = shufflelabel(DB);
figure, imageschs(dem,DB);
%hold on
%%
fprintf('saving temporary data for later use...')
DB_shuffled   = shufflelabel(DB);
save('Basin_shuffled_150m.mat','DB_shuffled','mask','x','y');

load Basin_shuffled_150m.mat
%%
%merge small basins with big basins from a minimum pixel number 3000 to 12000  
Basins = DB_shuffled.Z;
mask_fjords = mask; mask_fjords(mask==2) = 1; mask_fjords(mask~=1) = 0;
fprintf('merging small basins...')
for i = 3000:3000:12000
    mini_area = DB_shuffled.cellsize^2 * i;
    Basins = mergeBasins(Basins, DB_shuffled, mini_area);
    save(['Basin_mergered_',num2str(i),'.mat'], 'Basins');
%     figure,imagesc(mask_fjords, 'AlphaData',mask_fjords~=0.0),colormap(bone);
%     hold on;
%     imagesc(Basins,'AlphaData',Basins~=0.0),colormap(jet);
%     title(['mini = ',num2str(i)]);
end
%}
load Basin_shuffled_150m.mat
Basins_9000=load('Basin_mergered_9000.mat'); % minimum area of 135km^2

% get rid of the isolated ice patches
Basins=Basins_9000.Basins;
Basin_mask = bwareaopen(Basins,428494);
Basins(Basin_mask==0.0) = 0.0;
%Basins = Basins/100.0;

Basins(1575:1667,6169:6315)=0.0;
Basins(2375:2803,7227:7591)=0.0;
Basins(2874:3048,7625:7805)=0.0;
%NO
Basins(4178:4195,1104:1114)=0.0;
Basins(2243:2414,2426:2536)=0.0;

% redefine some basins
% now it is all manually 
% in the NO
Basins=GenMask_Basins_NO_150m(Basins);
Basins=GenMask_Basins_NE_150m(Basins);
Basins=GenMask_Basins_CE_150m(Basins);
Basins=GenMask_Basins_SE_150m(Basins);
Basins=GenMask_Basins_SW_150m(Basins);
Basins=GenMask_Basins_CW_150m(Basins);
Basins=GenMask_Basins_NW_150m(Basins);

save('Basins_150m.mat','Basins');
%%

fprintf('plotting...');
figure,imagesc(mask, 'AlphaData',(mask==1 | mask==2)),colormap(bone);
    hold on;
    imagesc(Basins,'AlphaData',Basins~=0.0),colormap(jet);
%{
figure, imagesc(I),colormap(jet);
figure, imagesc(xpsn),colormap(jet);
figure, imagesc(ypsn),colormap(jet);
figure, imagesc(mask),colormap(jet);

filename={'Basins_3000.Basins','Basins_6000.Basins','Basins_9000.Basins','Basins_12000.Basins'};
for i=1:4
    figure,imagesc(mask, 'AlphaData',(mask==1 | mask==2)),colormap(bone);
    hold on;
    tmp = eval(cell2mat(filename(i)));
    imagesc(tmp,'AlphaData',tmp~=0.0),colormap(jet);
    title(['mini = ',filename(i)]);
    clear tmp
end
%}