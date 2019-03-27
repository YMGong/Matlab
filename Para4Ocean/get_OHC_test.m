function [nearest_OHC, mean_OHC, fit_OHC,fit_err, nearest_xcoordi, nearest_ycoordi] = get_OHC_test(OHC, Ocean_coordix, Ocean_coordiy,...
                                                                    mouth_coordi, search_radius)
fprintf('This code does not search through the fjords!\n');
    %% create a bonding box and clip 
    r = search_radius;
    x = mouth_coordi(1);
    y = mouth_coordi(2);
    OHC_box = nan(size(OHC));
    step = abs(Ocean_coordix(2,1) - Ocean_coordix(1,1));
    while length(OHC_box(~isnan(OHC_box))) <= 100 % just to make sure there are at least 100 valide grid points 
        
        fprintf(['The search radius is ', num2str(r), '\n']);
        
        box = [(x-r/2), (x+r/2); (y-r/2), (y+r/2)];
        
        [x1, y1] = find(Ocean_coordix > box(1,1) & Ocean_coordix < box(1,2), 1, 'first');
        [x2, y2] = find(Ocean_coordix > box(1,1) & Ocean_coordix < box(1,2), 1, 'last');
        ocean_boxx = Ocean_coordix(x1:x2, y1:y2);
        ocean_boxy = Ocean_coordiy(x1:x2, y1:y2);
        OHC_box = OHC(x1:x2, y1:y2);
    
        [x1, y1] = find(ocean_boxy > box(2,1) & ocean_boxy < box(2,2), 1, 'first');
        [x2, y2] = find(ocean_boxy > box(2,1) & ocean_boxy < box(2,2), 1, 'last');
        ocean_boxx = ocean_boxx(x1:x2, y1:y2);
        ocean_boxy = ocean_boxy(x1:x2, y1:y2);
        xaxi = ocean_boxx(:,1);
        yaxi = ocean_boxy(1,:);
        
        OHC_box = OHC_box(x1:x2, y1:y2);
        
        ocean_boxx(isnan(OHC_box)) = nan;
        ocean_boxy(isnan(OHC_box)) = nan;
        r = r + 5*step; % Increase the search radius.
%                       The increased search radius only used when the modified OHC_box
%                       does not have enough valide grid points.
%                      length(OHC_box(~isnan(OHC_box)))
    end
    fprintf(['The fitting and averaging is done for ', num2str(length(OHC_box(~isnan(OHC_box)))), ' ocean grid points. \n']);
     figure, imagesc(Ocean_coordiy(1,:),Ocean_coordix(:,1),OHC),title('ocean_boxx');
     set(gca,'YDir','normal');
     hold on; 
     plot([box(2,1) ,box(2,1) ,box(2,2) , box(2,2), box(2,1)],...
             [box(1,1) ,box(1,2) ,box(1,2) , box(1,1), box(1,1)],...
         'k-', 'LineWidth', 3);  
    
    %% search the nearest ocean cell with in the box
  
    distance = sqrt((ocean_boxx - x) .^ 2 + (ocean_boxy - y) .^ 2);
    %figure, imagesc(distance)
    [idx_minD, idy_minD] = find(distance == min(min(distance)));
    % fit in a 2nd order polymonial
    [p, S, mu] = polyfit(distance(~isnan(distance)), OHC_box(~isnan(OHC_box)),2);
    %OHC_box(~isnan(OHC_box))
    [fit_OHC, fit_err] = polyval(p,min(min(distance)), S, mu);
    % or mean
    mean_OHC = nanmean(nanmean(OHC_box));
    % or use the nearest number directly
    nearest_OHC = OHC_box(idx_minD, idy_minD);
    nearest_xcoordi = ocean_boxx(idx_minD, idy_minD);
    nearest_ycoordi = ocean_boxy(idx_minD, idy_minD);
    OHC_box(idx_minD, idy_minD) = 2;
    %% plot and check
    figure, 
    h = imagesc(yaxi,xaxi,OHC_box);hold on;
    set(gca,'YDir','normal');
    set(h, 'AlphaData',~isnan(OHC_box));
    plot(y, x, '.r', 'MarkerSize',20);
    %{
    figure,  
    pcolorpsn(ocean_boxy,ocean_boxx, OHC_box,'meridian',-39);hold on;
    plot(nearest_xcoordi, nearest_ycoordi, '.r', 'MarkerSize',10);
    axis tight                    % gets rid of white space
%     mapzoompsn(82.22, -32.9735,'mapwidth',[800 500],'ne'); % zoom to the north
    xlabel('easting (m)','FontSize',15);
    ylabel('northing (m)','FontSize',15);
    %caxis([0 2.1])
    
    figure,  
    pcolorpsn(Ocean_coordiy,Ocean_coordix, OHC,'meridian',-39);
    axis tight                    % gets rid of white space
%     mapzoompsn(82.22, -32.9735,'mapwidth',[800 500],'ne'); % zoom to the north
    xlabel('easting (m)','FontSize',15);
    ylabel('northing (m)','FontSize',15);
    %}
    %}
end