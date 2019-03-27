function [data, Coordi]=readgeoTiff(filename)
if isfile(filename)
    data = geotiffread(filename);
    %%
    fprintf([filename,' information : \n']);
    data_info = geotiffinfo(filename);
    % convert the column and row indices of pixels to the world coordinate 
    [c,r] = pixcenters(data_info);
    [x, y] = meshgrid(c,r);
    [data_lat,data_lon] = projinv(data_info, x,y);
    Coordi(:,:,1) = data_lat;
    Coordi(:,:,2) = data_lon;
else
    error('file does not exist!\n');
end