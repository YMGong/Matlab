%************************************************************************%
% Save geotiff to nc 
% and interpolate z on 150m BedMachine grid
% 08/02/2019
% yongmei.gong@vub.be
%************************************************************************%

%clear all;
close all;

addpath /Users/yongmei/Documents/Academy/Coding/Matlab/Matlab_git/Functions
addpath /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_dem1978/data/0-data/G150AERODEM/DEM

%%
fprintf('Data loading ...\n');

% 1978 dem
[dem, dem_Coordi]=readgeoTiff('aerodem_1978_utm25.tif'); dem = rot90(dem, -1);
%dem1978(dem1978==min(min(dem1978)))=nan; % get rid of the noDatas
lat = dem_Coordi(:,:,1);lat = rot90(lat, -1);
lon = dem_Coordi(:,:,2);lon = rot90(lon, -1);

%%
fprintf('Saving...\n');
%infilename='aerodem_1978_utm25.mat';
outfilename='aerodem_1978_utm25';
%save(infilename,'lat','lon','dem1978');
row = size(dem,1);
colume = size(dem,2);
delete([outfilename,'.nc'])
system(['grdconvert /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_dem1978/data/0-data/G150AERODEM/DEM/',...
    outfilename, '.tif ', outfilename, '.nc '])
%system(['cp /Users/yongmei/Documents/Academy/Project/PARAMOUR/data_dem1978/data/0-data/G150AERODEM/DEM/', outfilename,' .'])

nccreate([outfilename,'.nc'],'lon','Dimensions',{'x' row 'y' colume},'Format','netcdf4');
ncwrite([outfilename,'.nc'],'lon',lon);
ncwriteatt([outfilename,'.nc'],'lon','units','degree_east');
ncwriteatt([outfilename,'.nc'],'lon','_CoordinateAxisType', "Lon");

nccreate([outfilename,'.nc'],'lat','Dimensions',{'x' row 'y' colume},'Format','netcdf4');
ncwrite([outfilename,'.nc'],'lat',lat);
ncwriteatt([outfilename,'.nc'],'lat','units','degree_north');
ncwriteatt([outfilename,'.nc'],'lat','_CoordinateAxisType', "Lat");

ncwriteatt([outfilename,'.nc'],'z','coordinates', 'lat lon');

system(['ncks -C -O -x -v x,y ' outfilename,'.nc ', outfilename, '.nc'])
%system('sh cdo_remapping.sh') it takes a lot of time!!
ncdisp(outfilename)

