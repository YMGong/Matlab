%*************************************************************************%
% The function search for the closest fjords points to the glacier ice front
% on the fjords model grid. The algrithem use breath-first search to pass 
% through the open fjords.
%      yongmei.gong@vub.be
%      19/03/2019
%*************************************************************************%
function [data_nearest, lon_nearest, lat_nearest] = get_ocean_data(ocean_data, ocean_lon, ocean_lat,...
    mask_fjords,BOcean_flag, fjords_lon, fjords_lat,...
         mouth_coordi, search_radius, fjordspointsThrh)
     %% intepolate ocean grid onto fjord grid
    xs = min(min(fjords_lon));  xe = max(max(fjords_lon));
    ys = min(min(fjords_lat)); ye = max(max(fjords_lat));
    stepD = abs(fjords_lon(2,1) - fjords_lon(1,1));
    stepV = 1;
    
    mask_ocean = ocean_data; 
    tmp_array = reshape(ocean_lon, size(ocean_lon,1)*size(ocean_lon,2),1);
    tmp_array(:,2) = reshape(ocean_lat, size(ocean_lat,1)*size(ocean_lat,2),1);
    tmp_array(:,3) = reshape(mask_ocean, size(mask_ocean,1)*size(mask_ocean,2),1);
    [mask_ocean, ~, ~]  = interp2array(xs,xe,ys,ye,tmp_array,stepV,stepD,'nearest');
    mask_ocean = rot90(mask_ocean);
    mask_ocean = fliplr(mask_ocean);
    %mask_ocean(~isnan(mask_ocean)) = 1;
    
     figure, imagesc(mask_ocean), colorbar;
    [nearest_xcoordi, nearest_ycoordi] = get_nearest_coordi(mask_ocean,...
         mask_fjords,BOcean_flag, fjords_lon, fjords_lat,...
         mouth_coordi, search_radius, fjordspointsThrh);
     
    tmp_array(:,4) = sqrt((double(tmp_array(:,1)) - double(nearest_xcoordi)) .^2 + ...
                          (double(tmp_array(:,2)) - double(nearest_ycoordi)) .^2);
    % the ocean points are nan...                  
    tmp_array(isnan(tmp_array(:,3)), 4) = nan;                            
                      
    [tmp_array(:,4), I] = sort(tmp_array(:,4));
    tmp_array(:,1:3) = tmp_array(I,1:3);
    lon_nearest = tmp_array(1, 1);
    lat_nearest = tmp_array(1, 2);
    data_nearest = tmp_array(1, 3);
    
    %%plot and check
    figure, 
    h = pcolor(ocean_lon, ocean_lat, ocean_data);
    hold on;
    set(h, 'EdgeColor', 'none');
    plot(lon_nearest, lat_nearest, '.k');
end 
function [nearest_xcoordi, nearest_ycoordi] = get_nearest_coordi(mask_ocean,...
         mask_fjords,BOcean_flag, fjords_lon, fjords_lat,...
         mouth_coordi, search_radius, fjordspointsThrh)
    %% create a bonding box for clipping 
    step = abs(fjords_lon(2,1) - fjords_lon(1,1)); % get the resolution
    r = round((search_radius/step)/2);
    fjords_lon = int64(fjords_lon);
    fjords_lat = int64(fjords_lat);
    [x, y] = find((mouth_coordi(1) ==  fjords_lon) &...
                  (mouth_coordi(2) ==  fjords_lat));
    %% search for at least 50 valid ocean grid points through breadth-first search
    mask_ocean_valid = 0;
    while length(mask_ocean_valid) <= fjordspointsThrh % OceanpointsThrh = 50
        %% clipping
        fprintf(['The search radius is ', num2str(r), '\n']);
        
        box = [(x-r), (x+r); (y-r), (y+r)]; % the clipping box
        mask_ocean_box = mask_ocean(box(1,1):box(1,2), box(2,1):box(2,2));
        ocean_boxx = fjords_lon(box(1,1):box(1,2), box(2,1):box(2,2));
        %ocean_boxx(isnan(mask_ocean_box)) = nan;
        ocean_boxy = fjords_lat(box(1,1):box(1,2), box(2,1):box(2,2));
        %ocean_boxy(isnan(mask_ocean_box)) = nan;
        mask_fjordsbox = mask_fjords(box(1,1):box(1,2), box(2,1):box(2,2));
        
        %% search for the connected ocean points connected with the fjords on the glacier model 
        nx = size(mask_fjordsbox, 1);
        ny = size(mask_fjordsbox, 2);
        [G, ~, ~] = get_search_graph(mask_fjordsbox, BOcean_flag); 
        [x_box, y_box] = find((mouth_coordi(1) ==  ocean_boxx) &...
                              (mouth_coordi(2) ==  ocean_boxy));  
        mouth_ID = toID(nx, ny, x_box, y_box);
        res = bfsearch(G, mouth_ID); % breadth-first search
        %% get the nearest ocean points on the ocean model grid along the search path 
        flag = 1;
        for i = 1:length(res)
            [resx,resy] = fromID(nx, ny, res(i));
            if ~isnan(mask_ocean_box(resx,resy)) 
                mask_ocean_valid(flag) = mask_ocean_box(resx,resy); 
                x_valid(flag) = resx; 
                y_valid(flag) = resy; 
%                 distance(flag) = sqrt(single((ocean_boxx(resx, resy) - mouth_coordi(1)) .^ 2 +...
%                                   (ocean_boxy(resx, resy) - mouth_coordi(2)) .^ 2) );
                flag = flag + 1; 
            end
           
        end
       r = r + 1; % Increase the search radius.
    end
    fprintf(['There are ', num2str(length(mask_ocean_valid)), ' valid ocean grid points in total. \n']);

    nearest_xcoordi = ocean_boxx(x_valid(1), y_valid(1));
    nearest_ycoordi = ocean_boxy(x_valid(1), y_valid(1));

    %% plot and check  
    figure, imagesc(mask_fjords),title('box'); hold on; 
    plot([box(2,1) ,box(2,1) ,box(2,2) , box(2,2), box(2,1)],...
         [box(1,1) ,box(1,2) ,box(1,2) , box(1,1), box(1,1)],...
             'k-', 'LineWidth', 3);  
     figure, imagesc(mask_ocean),title('box'); hold on; 
      plot([box(2,1) ,box(2,1) ,box(2,2) , box(2,2), box(2,1)],...
      [box(1,1) ,box(1,2) ,box(1,2) , box(1,1), box(1,1)],...
             'k-', 'LineWidth', 3);      
    figure, 
    h = imagesc(mask_ocean_box);hold on;
    set(h, 'AlphaData',~isnan(mask_fjordsbox));
    plot(x_box, y_box, '.r', 'MarkerSize',20);
    hold on;
    plot(y_valid(1), x_valid(1), '.k','MarkerSize',20); % x and y switched... stupid! 
end
%%
function [G, from, to] = get_search_graph(mask_search, GOcean_flag)
% the relative coordinate difference between the current and the surrounding cells
x_dif1 = [0 1 0 -1];
y_dif1 = [1 0 -1 0];
curIndex = 1;
nx = size(mask_search, 1);
ny = size(mask_search, 2);
% only search for the connection of the points indicated with GOcean_flag
   for x=2:nx-1
      for y=2:ny-1
           if mask_search(x,y) == GOcean_flag
               for i=1:4
                   if mask_search(x+x_dif1(i),y+y_dif1(i))== GOcean_flag
                       from(curIndex) = toID(nx,ny,x,y);
                       to(curIndex) = toID(nx,ny,x+x_dif1(i),y+y_dif1(i));
                       curIndex = curIndex+1;
                   end
               end
           end
      end
   end

   G = graph(from,to);                                                   
end
% convert indices for bfsearch
%%
function id = toID(nx, ny, x, y)
   if nx<ny
       id = (y-1)*ny + x;
   else
       id = (x-1)*nx + y;
   end
end
%%
function [x,y] = fromID(nx, ny, id)
   if nx<ny
       y = floor(id/ny)+1;
       x = id-(y-1)*ny;
       if x==0
           y = y-1;
           x = x+ny;
       end
   else
       x = floor(id/nx)+1;
       y = id-(x-1)*nx;
       if y==0
           x = x-1;
           y = y+nx;
       end
   end
end