%************************************************************************%
% Compare the masks
% C-GlORS data needs to be projected from regular lat-lon grid to 
% polar stereography grid and interpolated onto 1 km using cdo
% $cdo remapbil,$gridfile $infile $outfile
% 18/01/2019
% yongmei.gong@vub.be
%************************************************************************%
clear all;
%close all;
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_C-GLORS/C-GLORSV5/OHC/
%%
fprintf('Data loading ...\n');
config_greenland
% BedMachine mask
lonB = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lon'); lonB = rot90(lonB); % the x coordinates
latB = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lat'); latB = rot90(latB); % the y coordinates
[latBN,lonBN] = polarstereo_fwd(latB,lonB,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
% maskB = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','ima'); maskB = rot90(maskB); % 0= ocean; 1 = ice
% maskB(maskB~=0)=1;
maskB_fjords = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','cma'); maskB_fjords = rot90(maskB_fjords);
maskB_fjords(maskB_fjords==1)=2;
% C-GLOPRS NEMO mask
%{
maskN = ncread('C-GLORSv5_land_sea_mask_PolarStere.nc','lsm'); maskN=rot90(round(maskN)); %maskN = rot90(maskN); maskN=fliplr(maskN);
lonN = ncread('C-GLORSv5_land_sea_mask_PolarStere.nc','lon');lonN = rot90(lonN); 
latN = ncread('C-GLORSv5_land_sea_mask_PolarStere.nc','lat');latN = rot90(latN);
[latNN,lonNN] = polarstereo_fwd(latN,lonN,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
%}
maskN = ncread('C-GLORSv5_land_sea_mask.nc','lsm'); maskN=rot90(round(maskN)); %maskN = rot90(maskN); maskN=fliplr(maskN);
lonN = ncread('C-GLORSv5_land_sea_mask.nc','lon');lonN = rot90(lonN); 
lonN = wrapTo180(lonN);
latN = ncread('C-GLORSv5_land_sea_mask.nc','lat');latN = flipud(latN);
[lonN, latN] = meshgrid(lonN, latN);
%% get OHC
fprintf('Parameterization ...\n');
%mouth_coordi = [-1043000, -371000];
% for clipping 
xw = min(min(lonB)); xe = 0;%max(max(lonB));
ys = min(min(latB)); yn = max(max(latB));

mask_tmp = maskN;
mask_tmp((lonN<xw | lonN>xe) |...
    (latN<ys | latN>yn) ) = nan;
[is ,js] = find(~isnan(mask_tmp), 1, 'first');
[ie ,je] = find(~isnan(mask_tmp), 1, 'last');
maskN_gr =  maskN(is:ie, js:je);
lonN_gr = lonN(is:ie, js:je);
latN_gr = latN(is:ie, js:je);
[latN_gr,lonN_gr] = polarstereo_fwd(latN_gr,lonN_gr,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 

% mask_tmp(latN>ys & latN<yn) = 1;

%
mouth_coordi = [-834000 , 71000]; % lon; latn should use the first ocean grid point in front of the ice front on GISM 
search_radius = 50000; % the resl of nemo is about 1 km
fjordspointsThrh = 20; % only take the first 50 points
BOcean_flag = 0;
fprintf('get the data from the closet ocean points on the ocean model grid...\n');
data_tmp = maskN_gr; % should use the real OHC data 
%data_tmp(data_tmp == 0) = nan;
data_tmp(data_tmp == 0) = nan;
% we can choose if we want to use the mean, the nearest, or the fitted
% value (2nd order polymonial)
[data_nearest, lon_nearest, lat_nearest] = get_ocean_data(data_tmp, lonN_gr, latN_gr,...
                                             maskB_fjords,BOcean_flag, lonBN,latBN,...
                                             mouth_coordi, search_radius, fjordspointsThrh);

% [nearest_OHC, mean_OHC, fit_OHC,fit_err,...
%  nearest_xcoordi, nearest_ycoordi] = get_OHC_search(data_tmp,...
%  maskB_fjords,BOcean_flag, lonNN, latNN,...
%  mouth_coordi, search_radius, OceanpointsThrh);
%}
%%
%{
fprintf('Plotting ...\n');
%{
maskdiff=maskN;
maskdiff(maskB_fjords==1)=2;

% Plotting using the Arctic Mapping Tools
figure,  
%pcolorpsn(latN,lonN,maskdiff,'meridian',-39);
pcolorpsn(latN,lonN,maskN,'meridian',-39);hold on;
h = pcolorpsn(latN,lonN,maskB_fjords,'meridian',-39);
set(h, 'AlphaData',maskB_fjords==2,'FaceAlpha',.3);
caxis([0 2.1])
hold on;
contourpsn(latN,lonN, maskB_fjords,[1,1],'k','meridian',-39); 
colormap(jet);
axis tight                    % gets rid of white space
%greenland('k','meridian',-39) % plots black grounding line
xlabel('easting (m)','FontSize',15);
ylabel('northing (m)','FontSize',15);
%}
%}
