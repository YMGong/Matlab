%*************************************************************************%
% The function search for the closest ocean points to the glacier ice front
% on the ocean model grid. The algrithem use breath-first search to pass 
% through the open fjords.
%      yongmei.gong@vub.be
%      05/03/2019
%*************************************************************************%
 
function [nearest_OHC, mean_OHC, fit_OHC,fit_err,...
         nearest_xcoordi, nearest_ycoordi] = get_OHC_search(OHC,...
         mask_fjords,BOcean_flag, Ocean_coordix, Ocean_coordiy,...
         mouth_coordi, search_radius, OceanpointsThrh)
    %% create a bonding box for clipping 
    step = abs(Ocean_coordix(2,1) - Ocean_coordix(1,1)); % get the resolution
    r = round((search_radius/step)/2);
    Ocean_coordix = int64(Ocean_coordix);
    Ocean_coordiy = int64(Ocean_coordiy);
    [x, y] = find((mouth_coordi(1) ==  Ocean_coordix) &...
                  (mouth_coordi(2) ==  Ocean_coordiy));
    %% search for at least 50 valid ocean grid points through breadth-first search
    OHC_valid = 0;
    while length(OHC_valid) <= OceanpointsThrh % OceanpointsThrh = 50
        %% clipping
        fprintf(['The search radius is ', num2str(r), '\n']);
        
        box = [(x-r), (x+r); (y-r), (y+r)]; % the clipping box
        OHC_box = OHC(box(1,1):box(1,2), box(2,1):box(2,2));
        ocean_boxx = Ocean_coordix(box(1,1):box(1,2), box(2,1):box(2,2));
        %ocean_boxx(isnan(OHC_box)) = nan;
        ocean_boxy = Ocean_coordiy(box(1,1):box(1,2), box(2,1):box(2,2));
        %ocean_boxy(isnan(OHC_box)) = nan;
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
            if ~isnan(OHC_box(resx,resy)) 
                OHC_valid(flag) = OHC_box(resx,resy); 
                x_valid(flag) = resx; 
                y_valid(flag) = resy; 
                distance(flag) = sqrt(single((ocean_boxx(resx, resy) - mouth_coordi(1)) .^ 2 +...
                                  (ocean_boxy(resx, resy) - mouth_coordi(2)) .^ 2) );
                flag = flag + 1; 
            end
           
        end
       r = r + 1; % Increase the search radius.
    end
    %% get the fitted, mean or closest resultes 
    fprintf(['There are ', num2str(length(OHC_valid)), ' valid ocean grid points in total. \n']);
    fprintf(['But the fitting and averaging is only done for the first ', num2str(OceanpointsThrh), ' points. \n']);
    % fit in a 2nd order polymonial
    [p, S, mu] = polyfit(distance(~isnan(distance)), OHC_valid,2);
    [fit_OHC, fit_err] = polyval(p,distance(1), S, mu);
    % or mean
    mean_OHC = mean(mean(OHC_valid(1:OceanpointsThrh)));
    % or use the nearest single grid directly
    nearest_OHC = OHC_valid(1);
    nearest_xcoordi = ocean_boxx(x_valid(1), y_valid(1));
    nearest_ycoordi = ocean_boxy(x_valid(1), y_valid(1));
    OHC_box(x_valid(1), y_valid(1)) = 2;
    %% plot and check  
    figure, imagesc(OHC),title('box'); hold on; 
    plot([box(2,1) ,box(2,1) ,box(2,2) , box(2,2), box(2,1)],...
         [box(1,1) ,box(1,2) ,box(1,2) , box(1,1), box(1,1)],...
             'k-', 'LineWidth', 3);  
    figure, 
    h = imagesc(mask_fjordsbox);hold on;
    set(h, 'AlphaData',~isnan(mask_fjordsbox));
    plot(x_box, y_box, '.r', 'MarkerSize',20);
    hold on;
    plot(y_valid, x_valid, '.k'); % x and y switched... stupid! 
   
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