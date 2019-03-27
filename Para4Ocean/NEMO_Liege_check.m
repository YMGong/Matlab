%************************************************************************%
% Compare the data from NEMO
% 18/01/2019
% yongmei.gong@vub.be
%************************************************************************%
close all;
clear all;

addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_NEMO_Liege
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_BedMachine/Morlighem17_Datasets/New_dataset/
addpath /Users/yongmei/Library/'Application Support'/MathWorks/'MATLAB Add-Ons'/Toolboxes/'Arctic Mapping Tools'/'Arctic Mapping Tools'
%%
fprintf('data loading...')
bath = ncread('Bathy_ORCA24_Atl_Arctic_final_fusionne_modif_new_dom_modif_final_close_tide.nc', 'Bathymetry');
nav_lon = ncread('Bathy_ORCA24_Atl_Arctic_final_fusionne_modif_new_dom_modif_final_close_tide.nc', 'nav_lon');
nav_lat = ncread('Bathy_ORCA24_Atl_Arctic_final_fusionne_modif_new_dom_modif_final_close_tide.nc', 'nav_lat');

lonB = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lon'); lonB = rot90(lonB); % the x coordinates
latB = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lat'); latB = rot90(latB); % the y coordinates
maskB_fjords = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','cma'); maskB_fjords = rot90(maskB_fjords);
maskB_fjords(maskB_fjords==1)=2;
figure,  
%pcolorpsn(latN,lonN,maskdiff,'meridian',-39);
pcolorpsn(nav_lat,nav_lon,bath,'meridian',-39);hold on;
hold on;
contourpsn(latB,lonB, maskB_fjords,[1,1],'r','meridian',-39); 
colormap(jet);
axis tight                    % gets rid of white space
%greenland('k','meridian',-39) % plots black grounding line
xlabel('easting (m)','FontSize',15);
ylabel('northing (m)','FontSize',15);
mapzoompsn(82.22, -32.9735,'mapwidth',[800 500],'ne');