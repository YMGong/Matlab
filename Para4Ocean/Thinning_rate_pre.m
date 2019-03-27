%************************************************************************%
% Prepare for calculating thinning rate
% 08/02/2019
% yongmei.gong@vub.be
%************************************************************************%
clear all;
close all;

addpath /Users/yongmei/Library/'Application Support'/MathWorks/'MATLAB Add-Ons'/Toolboxes/'Arctic Mapping Tools'/'Arctic Mapping Tools'
addpath /Users/yongmei/Library/'Application Support'/MathWorks/'MATLAB Add-Ons'/Toolboxes/'crameri perceptually uniform scientific colormaps'/crameri
addpath /Users/yongmei/Documents/Academy/Coding/Matlab/Matlab_git/Functions
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_dem1978/data/0-data/G150AERODEM/DEM
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_GIMP/pub/DATASETS/nsidc0715_MEASURES_gimp_dem/reg
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_BedMachine/Morlighem17_Datasets/New_dataset/
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_BedMachine/Morlighem17_Datasets/Original_dataset/
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_MAR

%%
fprintf('parameterization...\n')
%load the constants for greenland
config_greenland    
yrS = [1979, 1981, 1985, 1987];
yrE = 2010;
rsl = 150;
%{
% for looping through old dems of the same year
year_olddem = {'1978','1981','1985','1987'};
path_olddem = '/Users/yongmei/Documents/Academy/Project/PARAMOUR/data_dem1978/data/0-data/G150AERODEM/DEM';

path_gimpdem = '/Users/yongmei/Documents/Academy/Project/PARAMOUR/data_GIMP/pub/DATASETS/nsidc0715_MEASURES_gimp_dem/reg'; 
sec_gimpdem = {'5_reg_30m_dem','4_reg_30m_dem','3_reg_30m_dem','2_reg_30m_dem','1_reg_30m_dem','0_reg_30m_dem'};
%}
%%
fprintf('Data loading ...\n');
load 1km_mask.mat % now use the mask with forjds openned
% 1=ice sheet, 3=bedrock, 4=ocean
mask = rot90(mask);
maskB_fjords = mask;
maskB_fjords(maskB_fjords==4)=nan;
maskB_ice = mask;
maskB_ice(maskB_ice==4 | maskB_ice==3)=nan;
% bedmechine mask
% maskB_fjords = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','cma'); maskB_fjords = rot90(maskB_fjords);
% maskB_ice = ncread('BedMachineGreenland_epsg3413_1kmBig_updatedGISM_Consistency_10mthk.nc','ima'); maskB_ice = rot90(maskB_ice); 
lonB = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lon'); lonB = rot90(lonB); % the x coordinates
latB = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lat'); latB = rot90(latB); % the y coordinates
[xB,yB] = polarstereo_fwd(latB,lonB,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
% convert lat/lon data for a polar stereographic system to decimal degrees
% on maps with negative numbers (-) for S and W. 
%}
%{
mask = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','mask'); mask = rot90(mask); 
% 0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land
maskB_fjords = mask; maskB_fjords(mask==3 | mask==4) = 0.0; 
xpsn = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','lon'); xpsn = rot90(xpsn); % the x coordinates
ypsn = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','lat'); ypsn = rot90(ypsn); % the x coordinates
[latB,lonB] = polarsereo_fwd(ypsn,xpsn,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
deltLat=latB(1,2)-latB(1,1);deltLon=lonB(1,1)-lonB(2,1);

% hard coded!!
% extend the area to include all the satellite images. Only work when the
% coordinates are inline
[latB_extend,lonB_extend] = meshgrid(latB(1,1)-(500*deltLat):deltLat:latB(end,end)+(500*deltLat),...
                                    lonB(1,1)+(500*deltLon):-deltLon:lonB(end,end)-(500*deltLon)); 
maskB_fjords_extend = NaN(size(latB_extend));
maskB_fjords_extend (501:end-500,501:end-500)= maskB_fjords;
%}

    % GIMP dem
%     [dem2011, dem2011_Coordi]=readgeoTiff('tile_3_5_reg_30m_dem.tif');
%     dem2011(dem2011==min(min(dem2011)))=nan; % get rid of the noDatas
%     dem2011(dem2011<=0)=nan;
    %figure, imagesc(dem2011),colorbar;
    %figure, imagesc(dem2011_Coordi(:,:,1)),colorbar;
%
%{
% 1979-1988 mean
smb80 = ncread('MAR80sSMBb.nc','smb'); smb80=rot90(smb80);
lon80 = ncread('MAR80sSMBb.nc','lon'); lon80=rot90(lon80);
lat80 = ncread('MAR80sSMBb.nc','lat'); lat80=rot90(lat80);
[x80,y80] = polarstereo_fwd(lat80,lon80,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
% smb80 = NaN(size(mask)); 
% smb80(int64(yB(1)-y80(1))/1000+1:end,1:size(x80,2))=smb80_tmp; %enlarge the domain to match Bedmachine
% clear smb80_tmp;
smb80=(rhoi/rhow)*smb80/1000.0; % convert to ice eq
smb80(maskB_ice~=1)=nan;
%}    
lonMAR = ncread('MARv3.9-yearly-ERA-Interim-1979.nc','LON'); lonMAR=rot90(lonMAR);
latMAR = ncread('MARv3.9-yearly-ERA-Interim-1979.nc','LAT'); latMAR=rot90(latMAR);
[xS,yS] = polarstereo_fwd(latMAR,lonMAR,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
smbMean = zeros(size(mask));

% interpolate to 150 m 
xs = -727925;  xe = 954625;
ys = -3459425; ye = -557675;
stepD = rsl;
stepV = 1;

for yr = 1979:1989
    smb_tmp = ncread(['MARv3.9-yearly-ERA-Interim-',num2str(yr),'.nc'],'SMB'); smb_tmp=rot90(smb_tmp);
   
    smb = NaN(size(mask)); 
    smb(int64(yB(1)-yS(1))/1000+1:end,1:size(xS,2))=smb_tmp; %enlarge the domain to match Bedmachine

    smb=(rhoi/rhow)*smb/1000.0; % convert to ice eq
    smb(maskB_ice~=1)=nan;
    smbMean = smbMean + smb;
end
smbMean = smbMean ./ length(1979:1989);
figure, h=imagesc(smbMean);colorbar,colormap(jet);set(h, 'AlphaData',~isnan(maskB_fjords));
title('mean');
clear smb_tmp smb
for i = 1%2:length(yrS) %1:length(yrS) 
%     %mar data
%     lonS = ncread(['MARv3.9-yearly-ERA-Interim-',num2str(yrS(i)),'.nc'],'LON'); lonS=rot90(lonS);
%     latS = ncread(['MARv3.9-yearly-ERA-Interim-',num2str(yrS(i)),'.nc'],'LAT'); latS=rot90(latS);
%     smb_tmp = ncread(['MARv3.9-yearly-ERA-Interim-',num2str(yrS(i)),'.nc'],'SMB'); smb_tmp=rot90(smb_tmp);

    %
    %calculating the anomalies 

%     [xS,yS] = polarstereo_fwd(latS,lonS,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
%     smbS = NaN(size(mask)); 
%     smbS(int64(yB(1)-yS(1))/1000+1:end,1:size(xS,2))=smb_tmp; %enlarge the domain to match Bedmachine
%     clear smb_tmp;
% 
%     smbS=(rhoi/rhow)*smbS/1000.0; % convert to ice eq
%     smbS(maskB_ice~=1)=nan;
    %
    %Deltasmb = zeros(size(smb80));
    Deltasmb = zeros(size(mask));
    for yr = yrS(i): yrE
        smb_tmp = ncread(['MARv3.9-yearly-ERA-Interim-',num2str(yr),'.nc'],'SMB'); smb_tmp=rot90(smb_tmp);
        smb = NaN(size(mask)); 
        
        smb(int64(yB(1)-yS(1))/1000+1:end,1:size(xS,2))=smb_tmp; %enlarge the domain to match Bedmachine
        smb=(rhoi/rhow)*smb/1000.0; % convert to ice eq m/yr
        smb(maskB_ice~=1)=nan;
        %Deltasmb = (smb - smb80) + Deltasmb;
        %Deltasmb = (smb - smbS) + Deltasmb;
        Deltasmb = (smb - smbMean) + Deltasmb;
    end
    fprintf('Interpolating...\n');
    if rsl == 150
        tmp_array = reshape(xB, size(xB,1)*size(xB,2),1);
        tmp_array(:,2) = reshape(yB, size(yB,1)*size(yB,2),1);
        tmp_array(:,3) = reshape(Deltasmb, size(Deltasmb,1)*size(Deltasmb,2),1);
        [Deltasmb, ~, ~]  = interp2array(xs,xe,ys,ye,tmp_array,stepV,stepD,'nearest');
    end
    figure, h=imagesc(Deltasmb);colorbar,colormap(jet);set(h, 'AlphaData',~isnan(Deltasmb));
    title(num2str(yrS(i)));
    save(['smb_anom_Mean_',num2str(yrS(i)),'_',num2str(rsl),'.mat'],'Deltasmb');
end
%{    
for i=2:length(year_olddem)
    
%     files = dir (fullfile(path_olddem,['*',cell2mat(year_olddem(i)),'*.tif']));
%     L = length (files);
%     for k=1
%         filename=strcat(path_olddem,'/',files(k).name);
        [newmap_lon, newmap_lat, demold_gr] = patching_dems(path_olddem,cell2mat(year_olddem(i)), lonB_extend, latB_extend);
        save(['olddem_',cell2mat(year_olddem(i)),'.mat'], 'demold_gr');
    %     fprintf('Plotting...\n');
        figure, h=imagesc(demold_gr);
                colormap(jet);hold on;
                set(h, 'AlphaData',~isnan(demold_gr));
    %           contour(maskB_fjords_extend, 1);
%     end
end

for j=1:length(sec_gimpdem)

    [newmapgimp_lon, newmapgimp_lat,gimpdem_gr] = patching_dems(path_gimpdem, cell2mat(sec_gimpdem(j)), lonB_extend, latB_extend);
    save(['gimpdem_',cell2mat(sec_gimpdem(j)),'.mat'], 'gimpdem_gr');
%     fprintf('Plotting...\n');
    figure, h=imagesc(gimpdem_gr);
             colormap(jet);hold on;
             set(h, 'AlphaData',~isnan(gimpdem_gr));
           contour(maskB_fjords_extend, 1);
end
%}
%{
%%
fprintf('Plotting...\n');
   figure,
   h = imagesc(x80(1,:), y80(:,1), smb80);
   set(h, 'AlphaData',maskB_fjords==1);set (gca,'Ydir','normal');
   %axis([-9e5 1e6 -3.4e6 -0.6e6])
   colorbar,colormap(jet);%lcolorbar('m i.e./yr')

    
     figure, h = pcolorpsn(dem1978_Coordi(:,:,1), dem1978_Coordi(:,:,2), dem1978);colorbar; % check with the Arctic Mapping Tools
     set(h, 'AlphaData',~isnan(dem1978));
     crameri('devon');shadem(1);  
    contourpsn(latB,lonB, maskB_fjords,[1,1],'k'); 
    axis tight                    % gets rid of white space
    %greenland('k','meridian',-39) % the build in greenland contour is off! DO
    %NOT USE
    xlabel('easting (m)','FontSize',15);
    ylabel('northing (m)','FontSize',15);
    mapzoompsn(82.22, -32.9735,'mapwidth',[800 500],'ne'); % zoom to the north
    %mapzoompsn(82.22, -32.9735,'mapwidth',[1800 500],'ne'); % zoom to the entire north
    saveas(gcf,'dem1978_NO_1stpatch.tif','tif')

    figure, h = pcolorpsn(dem2011_Coordi(:,:,1), dem2011_Coordi(:,:,2), dem2011); colorbar;% check with the Arctic Mapping Tools
     set(h, 'AlphaData',~isnan(dem2011));
     crameri('devon');shadem(1); 
    contourpsn(latB,lonB, maskB_fjords,[1,1],'k'); 
    axis tight                    % gets rid of white space
    %greenland('k','meridian',-39) % the build in greenland contour is off! DO
    %NOT USE
    xlabel('easting (m)','FontSize',15);
    ylabel('northing (m)','FontSize',15);
    mapzoompsn(82.22, -32.9735,'mapwidth',[800 500],'ne'); % zoom to the north
    %mapzoompsn(82.22, -32.9735,'mapwidth',[1800 500],'ne'); % zoom to the entire north
    saveas(gcf,'dem2011_NO_3_5.tif.tif','tif')
    %}

