function [basemap_lon, basemap_lat, patching_basemap] = patching_dems_test(filepath, keywords, basemap_lon, basemap_lat)
    %load the constants for greenland
    constants_greenland
    % get the file names that has the keywords
    files = dir (fullfile(filepath,['*',keywords,'*.tif']));
    L = length (files);
   %creat a big greenland map for patching data
        patching_basemap = zeros(size(basemap_lon));
        patching_basemap(patching_basemap==0)=NaN;
    % for interpolation
    stepD = basemap_lon(1,1)-basemap_lon(2,1); % stepping through the regular array for nearest
                     % neighbour calculations (in practice this gives grid
                     % size in metres, you probably want somewhere between
                     % 100m - 500m for Austfonna)
    stepV = 1; % stepping through the scattered data (set this to 1 to
                   % use all the data, or a higher integer to skip some
                   % data, which will speed up this script at the cost of
                   % accuracy)
    size_original = size(basemap_lon);
    for i=1:L
        %fprintf(files(i).name);
        [data, data_Coordi]=readgeoTiff(strcat(filepath,'/',files(i).name));
        %dem1978(dem1978==min(min(dem1978)))=nan; % get rid of the noDatas
        data(data<=0)=nan;
        %%
        fprintf('aligning data together...\n');
        
        [lat_data, lon_data] = polarstereo_fwd(data_Coordi(:,:,1),data_Coordi(:,:,2),EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
%         figure, imagesc(lat_data);colormap(jet);colorbar;
%         figure, imagesc(lon_data);colormap(jet);colorbar;
        
        %lat_olddem = single(lat_olddem); lon_olddem = single(lon_olddem);
        % for interp2array
        % get the closest corner points of dem1978 on demB
        
        close_point = find(basemap_lat<=min(min(lat_data)));
        % check if the patch is within the basemap, otherwise increase the
        % size of the basemap
        while isempty(close_point)
            basemap_lat = [basemap_lat(:,1)-stepD, basemap_lat];
            close_point = find(basemap_lat<=min(min(lat_data)));
        end
        [~,ymin] = find(basemap_lat==max(basemap_lat(close_point)));
        ymin = unique(ymin);
        
        close_point = find(basemap_lon(:,ymin)>=max(max(lon_data)));
        while isempty(close_point) 
            basemap_lon = [basemap_lon(1,:)+stepD; basemap_lon];
            close_point = find(basemap_lon(:,ymin)>=max(max(lon_data)));
        end
        %min(basemap_lon(close_point))
        [xmin, ~] = find(basemap_lon==min(basemap_lon(close_point)));
        xmin = unique(xmin);

        close_point = find(basemap_lat>=max(max(lat_data)));
        while isempty(close_point)
            basemap_lat = [basemap_lat, basemap_lat(:,1)+stepD];
            close_point = find(basemap_lat>=max(max(lat_data)));
        end
        [~,ymax] = find(basemap_lat==min(basemap_lat(close_point)));
        ymax = unique(ymax);

        close_point = find(basemap_lon(:,ymax)<=min(min(lon_data)));
        while isempty(close_point)
            basemap_lon = [basemap_lon; basemap_lon(1,:)-stepD];
            close_point = find(basemap_lon(:,ymax)<=min(min(lon_data)));
        end
        [xmax, ~] = find(basemap_lon==max(basemap_lon(close_point)));
        xmax = unique(xmax);

        xs = basemap_lat(xmin,ymin);  xe = basemap_lat(xmax,ymax);
        ys = basemap_lon(xmax,ymax); ye = basemap_lon(xmin,ymin);

        tmp_array = reshape(lat_data, size(lat_data,1)*size(lat_data,2),1);
        tmp_array(:,2) = reshape(lon_data, size(lon_data,1)*size(lon_data,2),1);
        tmp_array(:,3) = reshape(data, size(data,1)*size(data,2),1);
        [demold_150m, ~, ~]  = interp2array(xs,xe,ys,ye,tmp_array,stepV,stepD,'linear');
        
        patching_basemap(xmin:xmax, ymin:ymax) = demold_150m; %patch it on to the big greenland icesheet grid

        clear tmp_array close_point 
    end
    if (size(basemap_lat,1)~=size_original(1) || size(basemap_lat,2)~=size_original(2))
        fprintf('*** The basemap has been changed ! ***\n');
    end
end