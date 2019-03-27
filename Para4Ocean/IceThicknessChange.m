%************************************************************************%
% process the GIMP DEM time series
% 14/01/2018
% yongmei.gong@vub.be
%************************************************************************%
clear all;
close all;
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_GIMP/pub/DATASETS/nsidc0715_MEASURES_gimp_dem/reg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Data loading ...\n');

[dem,R] = geotiffread('tile_0_5_reg_2012_2_30m_dem.tif');
figure,imagesc(dem)%,colormap(jet)
% mapshow(dem,R);
% axis image off