function [centerlinex, centerliney, centerline_distance, surface, bed, dynthinning] =...
    find_centerline(xcoor, ycoor, surfaceDEM, bedDEM,r_dyn,...
    interpDistThreshold, sample_distance, mask, step, terminus, PointE,...
    maxdelta, edge_cutoff, distanceThresh, varargin)
% find the centerline and interpolate data along center line    

switch varargin{1}
    case 'find_centerline_twopoints'
        centerPoint = find_centerline_twopoints(mask, xcoor, ycoor,...
        terminus, PointE, maxdelta, step);
    case 'find_centerline_onepoint'
        centerPoint = find_centerline_onepoint(mask,  xcoor, ycoor,...
        terminus, step, edge_cutoff, maxdelta);
    otherwise
    error(['Unexpected option: ' varargin]);    
end

    for i = 1:length(centerPoint)
        x(i) = xcoor(centerPoint(1, i), centerPoint(2, i));
        y(i) = ycoor(centerPoint(1, i), centerPoint(2, i));
        surface(i)= surfaceDEM(centerPoint(1, i), centerPoint(2, i));
        bed(i)= bedDEM(centerPoint(1, i), centerPoint(2, i));
        dynthinning(i) = r_dyn(centerPoint(1, i), centerPoint(2, i));
        if isnan(dynthinning(i))
            centerPoint(1, i)
            centerPoint(2, i)
        end
    end
    
    
    % calculate the accumulative distance
    d = ((x(2:end) - x(1:end-1)).^2 + (y(2:end) - y(1:end-1)).^2).^(1/2);
    d = [0, cumsum(d)];
    d = d(d < distanceThresh);
    x = x(d < distanceThresh);
    y = y(d < distanceThresh);
    surface = surface(d < distanceThresh);
    bed = bed(d < distanceThresh);
    dynthinning = dynthinning(d < distanceThresh);
    
    % resample 
    % extend the centerline to include some poi
    centerline_distance = -100* sample_distance : sample_distance :d(end);
    % for checking is the points are too far away from each other
    % we might not need it due to the way we sample the center line points
     if isnan(interpDistThreshold) 
        for i = 1:length(centerline_distance)
            distance(i) = min(abs(centerline_distance(i) - d)); 
        end
        centerline_distance(distance > interpDistThreshold) = nan;
     end
    % extrapolate to include some points that are in the ocean 
    % for getting OHC from the ocean model
    centerlinex = interp1(d,x,centerline_distance, 'linear', 'extrap');
    centerliney = interp1(d,y,centerline_distance,'linear', 'extrap');
    surface = interp1(d,surface,centerline_distance);
    bed = interp1(d,bed,centerline_distance);
    dynthinning = interp1(d,dynthinning,centerline_distance);
    % check the plots
        
%     figure, h=imagesc(mask);colormap(jet);colorbar;
%     set(h, 'AlphaData',~isnan(mask));
%     title('old')
%     hold on;plot(centerPoint(2,:), centerPoint(1,:),'.y');
    
    

    figure, h=imagesc(ycoor(1,:), xcoor(:,1), mask);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(mask));
    set(gca,'YDir','normal')
    title('Uninterpolated centerline points')
    hold on
    plot(y, x, '.y');
    set(gca,'FontSize',15)
    
    figure, h=imagesc(ycoor(1,:), xcoor(:,1), mask);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(mask));
    set(gca,'YDir','normal')
    title('Interpolated centerline points')
    hold on
    plot(centerliney, centerlinex, '.k');
    set(gca,'FontSize',15)
    
    figure,
    test_plot('Profile along centerline',centerline_distance(centerline_distance>=0),...
    'surface', surface(~isnan(surface)),...
    'bed', bed(~isnan(bed)));
    
    figure,
    plot(centerline_distance(centerline_distance>=0), dynthinning(centerline_distance>=0),...
        'b', 'LineWidth',2);
    title('Profile along centerline');
    legend('r_dyn');
    xlabel('Distance along certerline (m)');ylabel('dynamic thinning (m/yr)');
    set(gca,'FontSize',15)
%}
end
function centerPoint = find_centerline_onepoint(mask,  xcoor, ycoor,...
        terminus, step, edge_cutoff, maxdelta)
     if max(max(isnan(mask)))
         error('No data points need to be set to 0 !\n')
     end
     [rowi, coli] = find((terminus(1) ==  int64(xcoor)) &...
                  (terminus(2) ==  int64(ycoor)));
     row = rowi;
     col_len = size(mask,2);
     centerPoint = [rowi; coli];
    while max(mask(row+step,:))
        %max(bwareaopen(mask(3611,:),3))
        row = row+step;
        firstI = find(bwareaopen(mask(row,:),edge_cutoff)~=0,1,'first'); % get rid of the disconnected points in each row
        if ~isempty(firstI) 
            %lastI = find(mask(row,:)~=0,1,'last');
            lastI = firstI;
            while mask(row,lastI)
                lastI = lastI+1;
            end
            lastI = lastI-1;
            delta = round(abs(lastI - firstI)/2);
            if delta > maxdelta
                delta = maxdelta;
            end
            centerPoint = [centerPoint,...
                    [row;  firstI + delta]];
        end
    end
end
function centerPoint = ...
    find_centerline_twopoints(mask, xcoor, ycoor, terminus, PointE, maxdelta, step, varargin)
%% the starting point and ending point
% default
PLine_step = 1;
line_threshold = 100;
for i = 1:2:length(varargin)
            switch varargin{i}
                case 'PLine_step'
                      PLine_step = varargin{i+1};      
                case 'line_threshold'
                    line_threshold = varargin{i+1};
                otherwise
                      error(['Unexpected option: ' varargin{i}]);
            end
end
[i_s, j_s] = find((terminus(1) ==  int64(xcoor)) &...
                  (terminus(2) ==  int64(ycoor)));
              
[i_e, j_e] = find((PointE(1) ==  int64(xcoor)) &...
                  (PointE(2) ==  int64(ycoor)));              

pt1 = [i_s, j_s];
pt2 = [i_e, j_e];

n = max(j_e -j_s,line_threshold);
%n = 100;
t = linspace(0,1,n+1); 
%t = t(2:(end-1)); 
t = t(1:end); 
v = pt2 - pt1;
x = pt1(1) + t*v(1);    % p(t) = p1 + t*(p2-p1)
y = pt1(2) + t*v(2);
figure, imagesc(mask), hold on;
plot(y, x, 'r'); hold on;
% find the perpendicular lines
width = 1;
inum = size(mask, 1);
jnum = size(mask, 2);
centerline_r =  i_s;
centerline_c =  j_s;
for r= 4:step:n
    v = width*v / norm(v);
    
    r_ll = max(round(x(r)+v(2)),1); 
    r_ll = min(r_ll, inum); 
    
    c_ll = max(round(y(r)-v(1)), 1);
    c_ll = min(c_ll, jnum);
    
    r_ur = max(round(x(r)-v(2)), 1);
    r_ur = min(r_ur, inum);
    
    c_ur = max(round(y(r)+v(1)), 1);
    c_ur = min(c_ur, jnum);
    
    while (( mask( r_ll, c_ll ) ~=0 ) || ...
             ( mask( r_ur, c_ur ) ~=0 ) )
         
        width=width+1;
        v = width*v / norm(v);
        r_ll = max(round(x(r)+v(2)),1); 
        r_ll = min(r_ll, inum); 

        c_ll = max(round(y(r)-v(1)), 1);
        c_ll = min(c_ll, jnum);

        r_ur = max(round(x(r)-v(2)), 1);
        r_ur = min(r_ur, inum);

        c_ur = max(round(y(r)+v(1)), 1);
        c_ur = min(c_ur, jnum);
    
    end 
    line([y(r)-v(1), y(r)+v(1)], [x(r)+v(2), x(r)-v(2)]);
    pbaspect([1 1 1]);
    slope = (r_ur-r_ll) / (c_ur-c_ll); %scalar
    if c_ll <= c_ur
        c_new = c_ll:PLine_step:c_ur;
        r_new = round(slope * (0:PLine_step:(c_ur - c_ll)) + r_ll); 
    else
        c_new = c_ll:-PLine_step:c_ur;
        r_new = round(slope * (0:-PLine_step:(c_ur - c_ll)) + r_ll); 
    end
    %counter = 1;
    FirstI = 0;
    for k = 1:length(c_new)
        if mask(r_new(k), round(c_new(k))) ~=0 && FirstI == 0
            FirstI = 1;
            mask(r_new(k), round(c_new(k))) = 2;
            r_first = r_new(k);
            c_first = c_new(k);
            break;
            %counter = counter + 1;
        end
    end
    LastI = 0;
    
    for h = k:length(c_new)
        if mask(r_new(h), round(c_new(h))) == 0 && LastI == 0         
            LastI = 1;
            mask(r_new(h-1), round(c_new(h-1))) = 2;
            r_last = r_new(h-1); % the point before the 0 point 
            c_last = c_new(h-1);
            break;
            %counter = counter + 1;
        end
    end

    centerline_r = [centerline_r, r_first - min(maxdelta, round((r_first - r_last) / 2) )];
    centerline_c = [centerline_c, c_first + min(maxdelta, round((c_last - c_first) / 2) )];
%     [centerline_r, I] = sort(centerline_r);
%     centerline_c = centerline_c(I);
    centerPoint = [centerline_r; centerline_c];
    
%     mask(r_last + round((r_first - r_last) / 2),...
%         c_first + round((c_last - c_first) / 2)) = 3;

    clear r_tmp c_tmp
end
figure, imagesc(mask); pbaspect([1 1 1]);colorbar;
%set(h, 'alphadata', mask~=0);
hold on;
plot(y, x, 'r');
hold on;
plot(centerline_c, centerline_r, 'b.');
end