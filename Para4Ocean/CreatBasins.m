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
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_BedMachine/Morlighem17_Datasets/New_dataset/
addpath /Users/yongmei/Library/'Application Support'/MathWorks/'MATLAB Add-Ons'/Toolboxes/'Arctic Mapping Tools'/'Arctic Mapping Tools'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Data loading ...\n');

%ncdisp('sur.05_g.grd')
%I = ncread('sur.01.cdf','sur'); I = rot90(I); % the surface elevation
% mask = ncread('mask05.cdf','z'); % 0 = shelf; 1 = sheet; 3 = bare rock; 4 = ocean
% mask = ncread('ima.01.cdf','ima'); mask = rot90(mask); % 0= ocean; 1 = ice
% mask(mask~=0)=1;

%11/01/2018 below I use the data with minimum 10m ice thickness instead of
%50m
I = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','sur'); I = rot90(I); % the surface elevation
bed = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','bed'); bed = rot90(bed);
figure, imagesc(bed)
xpsn = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lon'); xpsn = rot90(xpsn); % the x coordinates
ypsn = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lat'); ypsn = rot90(ypsn); % the y coordinates
mask = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','ima'); mask = rot90(mask); % 0= ocean; 1 = ice
mask(mask~=0)=1;
mask_fjords = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','cma'); mask_fjords = rot90(mask_fjords);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Parameterization ...\n');
%load the constants for greenland
constants_greenland 
saveGridfile = 0; % 0 = do not save the GridObj file created by mat2gridobj
resl = 1000; %1km grid size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Find the basins ...\n');
[x,y] = polarstereo_fwd(ypsn,xpsn,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
% convert lat/lon data for a polar stereographic system to decimal degrees
% on maps with negative numbers (-) for S and W. 
xllcorner = x(end,1); % ideally x coordinate of the Center Point of the lower left Corner Grid Cell
yllcorner = y(end,1); % ideally y coordinate of the Center Point of the lower left Corner Grid Cell
%[xutm,yutm] = ll2utm(x,y,'wgs84',19);
dem = mat2gridobj(I,mask, xllcorner, yllcorner, resl, saveGridfile); %creat a GRIDobj
hillshade(dem)
FD  = FLOWobj(dem,'preprocess','c');
DB   = drainagebasins(FD);
%DB   = shufflelabel(DB);
%imageschs(dem,DB);
%hold on

%merge small basins with big basins from a minimum pixel size 500 to 2000 (area from 500km to 2000km) 
Basins = DB.Z;
for i = 500%:1000:2500
    mini_area = DB.cellsize^2 * i;
    Basins = mergeBasins(Basins, DB, mini_area);
%     figure,imagesc(mask_fjords, 'AlphaData',mask_fjords~=0.0),colormap(bone);
%     hold on;
%     imagesc(Basins,'AlphaData',Basins~=0.0),colormap(jet);
%     title(['mini = ',num2str(i)]);
end
% get rid of the isolated ice patches

Basins(227:351,958:1455)=nan;
Basins(335:417,1104:1218)=nan;
Basins(308:330,964:981)=nan;
Basins(353:365,932:973)=nan;
Basins(366:369,936:948)=nan;
Basins(366:371,956:965)=nan;
Basins(351:383,593:630)=nan;


Basins(326:340,718:738)=nan;
Basins(403:404,674:676)=nan;
Basins(383:401,582:601)=nan;
Basins(382:461,372:434)=nan;
Basins(422:436,459:468)=nan;
Basins(446:502,329:390)=nan;

Basins(742:762,138:155)=nan;
Basins(878:894,234:253)=nan;

Basins(1474:1476,495:498)=nan;
Basins(1653:1694,396:439)=nan;
Basins(1596,463)=nan;
Basins(2039:2104,382:485)=nan;
Basins(1566:1750,379:481)=nan;


Basins(2710,789)=nan;
Basins(2222:2231,920:943)=nan;
Basins(2222:2231,920:943)=nan;
Basins(2126:2137,1003:1011)=nan;
Basins(2054:2073,1030:1036)=nan;


Basins(495:539,1166:1221)=nan;
Basins(235:416,1060:1242)=nan;
Basins(302:371,954:979)=nan;

Basins(1328:1496,1334:1405)=nan;
Basins(1128:1146,1272:1286)=nan;
Basins(1113:1118,1298:1306)=nan;
Basins(1106:1116,1313:1320)=nan;
Basins(1405:1413,1262:1269)=nan;
Basins(2012:2023,1056:1066)=nan;
Basins(319:341,714:733)=nan;
Basins(369:373,939:962)=nan;
Basins(350:354,934:964)=nan;

Basins(485:524,1159:1186)=nan;
Basins(537:545,1171:1233)=nan;
Basins(539:558,1219:1234)=nan;

Basins(841:851,1221:1230)=nan;
Basins(1126,1339)=nan;
Basins(1131:1134,1231:1234)=nan;
Basins(1059:1064,1309:1315)=nan;
Basins(1063:1069,1273:1279)=nan;
Basins(1111:1113,1260:1265)=nan;

Basins(1426,1333)=nan;
Basins(1710:1732,1315:1352)=nan;
Basins(1758:1770,1305:1309)=nan;
Basins(1776:1786,1290:1293)=nan;
Basins(1763:1772,1273:1279)=nan;
Basins(1749:1766,1209:1222)=nan;
Basins(1878:1885,1180:1186)=nan;

Basins(1963:1969,1140:1150)=nan;
Basins(2303:2316,878:892)=nan;
Basins(2355:2360,877:883)=nan;
Basins(2481:2484,809:814)=nan;
Basins(2707:2739,775:793)=nan;
Basins(2603:2616,797:802)=nan;
Basins(1427,1333)=nan;

Basins(2304:2312,487:495)=nan;
Basins(2285:2288,498:501)=nan;
Basins(1482:1487,512:515)=nan;
Basins(2063:2107,372:450)=nan;
Basins(1472:1480,492:501)=nan;
Basins(1455:1458,501:504)=nan;

% redefine some basins
% now it is all manually 
% in the NO
Basins=GenMask_Basins_NO_1km(Basins);
% in the NE
Basins=GenMask_Basins_NE_1km(Basins);
%in the CE
Basins=GenMask_Basins_CE_1km(Basins);
%in the SE
Basins=GenMask_Basins_SE_1km(Basins);
%in the SW
Basins=GenMask_Basins_SW_1km(Basins);
%in the CW
Basins=GenMask_Basins_CW_1km(Basins);
%in the NW
Basins=GenMask_Basins_NW_1km(Basins);


%{

Basins(Basins==2310)=433;
Basins(Basins==1635|Basins==1454|Basins==491|Basins==377)=960;
Basins(Basins==2666)=782;
Basins(Basins==1635|Basins==1454|Basins==491|Basins==377)=960;

temp=Basins(349:390, 901:928);
temp(temp==1772) = 971;
Basins(349:390, 901:928)=temp;
temp=Basins(390:416,977:1002);
temp(temp==1772) = 1042;
Basins(390:416,977:1002)=temp;

Basins(396:404,976)=1042;
Basins(Basins==1772)=2362;

Basins(Basins==2298)=2414;

Basins(Basins==863)=1738;

temp=Basins(546:650,1179:1226);
temp(temp==2414) = 1738;
Basins(546:650,1179:1226)=temp;

temp=Basins(712:739,1190:1257);
temp(temp==2668) = 555;
Basins(712:739,1190:1257)=temp;

temp=Basins(373:750,1190:1215);
temp(temp==2668) = 555;
Basins(373:750,1190:1215)=temp;

Basins(Basins==1416|Basins==2668)=2755;

temp=Basins(1023:1091,1307:1390);
temp(temp==1411) = 511;
Basins(1023:1091,1307:1390)=temp;

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Plotting ...\n');
figure,imagesc(mask_fjords, 'AlphaData',mask_fjords~=0.0),colormap(bone);
hold on;
imagesc(Basins,'AlphaData',Basins~=0.0);colormap(jet);

% shuffling the color to make the basins more distinguishable 
DB.Z=Basins;
DB   = shufflelabel(DB); 
%mask_fjords(mask_fjords~=1)=nan;
DB.Z(Basins==0.0)=nan;
figure,  
%pcolorpsn(ypsn,xpsn,mask_fjords,'meridian',-39);
%hold on;
h = pcolorpsn(ypsn,xpsn,DB.Z,'meridian',-39);hold on;
%set(h,'FaceAlpha',  'texturemap', 'AlphaDataMapping', 'none', 'AlphaData',DB.Z~=0.0);
contourpsn(ypsn,xpsn, mask_fjords,[1,1],'k','meridian',-39); hold on
axis tight                    % gets rid of white space
%greenland('k','meridian',-39) % plots black grounding line
xlabel('easting (m)','FontSize',15);
ylabel('northing (m)','FontSize',15);

figure,imagesc(mask_fjords, 'AlphaData',mask_fjords~=0.0),colormap(bone);
hold on;
h=imagesc(DB.Z,'AlphaData',Basins~=0.0);colormap(jet);
% title(['mini = ',num2str(i)]);

bed(mask_fjords==1)=nan;
figure,
pcolorpsn(ypsn,xpsn,bed,'meridian',-39);hold on;
contourpsn(ypsn,xpsn, mask_fjords,[1,1],'k','meridian',-39); hold on
axis tight                    % gets rid of white space
%greenland('k','meridian',-39) % plots black grounding line
xlabel('easting (m)','FontSize',15);
ylabel('northing (m)','FontSize',15);
%NE
%set(h,'alphadata',Basins==243|Basins==326|Basins==2480|Basins==2914|Basins==368)
%{
fprintf('Merge basins according to Porter et al.,2018...\n');

Basins = D.Z;
Basins(mask==0)=nan;
Basins_N = Basins(200:1100, :);
Basins_N=GenMask_Basins_N_1km(Basins_N);
figure, imagesc(Basins_N, 'AlphaData',Basins_N~=0.0),colormap('jet');colorbar;

Basins_Thule = Basins(550:900,100:420);
Basins_Thule=GenMask_Basins_Thule_1km(Basins_Thule);
figure, imagesc(Basins_Thule, 'AlphaData',Basins_Thule~=0.0),colormap('jet');colorbar;

Basins_W = Basins(741:2707, 252:1014);
Basins_W=GenMask_Basins_W_1km(Basins_W);
figure, imagesc(Basins_W, 'AlphaData',Basins_W~=0.0),colormap('jet');colorbar;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Plotting ...\n');
figure, imagesc(Basins_N, 'AlphaData',Basins_N~=0),colormap('jet');colorbar;
%figure, imageschs(dem,D,'falsecolor',[1,1,1]);
%figure, h=imageschs(dem,D,'falsecolor',[1,1,1]);
%figure, imagesc(mask),colorbar;
%figure, imagesc(dem),colorbar;
%}