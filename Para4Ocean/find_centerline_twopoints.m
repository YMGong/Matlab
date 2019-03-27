function find_centerline_twopoints(mask, xcoor, ycoor, terminus, PointE, maxdelta, varargin)
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
                      error(['Unexpected option: ' varargin{i}])
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
%% find the perpendicular lines
width = 1;
inum = size(mask, 1);
jnum = size(mask, 2);
centerline_r =  i_s;
centerline_c =  j_s;
for r= 5:5:n
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
        if mask(r_new(h), round(c_new(h))) ==0 && LastI == 0
            
            LastI = 1;
            mask(r_new(h-1), round(c_new(h-1))) = 2;
            r_last = r_new(h-1); % the point before the 0 point 
            c_last = c_new(h-1);
            break;
            %counter = counter + 1;
        end
    end
    length(r_last +  round((r_first - r_last) / 2))
    length(c_first + round((c_last - c_first) / 2))
%     centerline_r = [centerline_r, r_first -  round((r_first - r_last) / 2) ];
%     centerline_c = [centerline_c, c_first + round((c_last - c_first) / 2) ];
    centerline_r = [centerline_r, r_first - min(maxdelta, round((r_first - r_last) / 2) )];
    centerline_c = [centerline_c, c_first + min(maxdelta, round((c_last - c_first) / 2) )];
    
%     mask(r_last + round((r_first - r_last) / 2),...
%         c_first + round((c_last - c_first) / 2)) = 3;

    clear r_tmp c_tmp
end
figure, h=imagesc(mask); pbaspect([1 1 1]);colorbar;
%set(h, 'alphadata', mask~=0);
hold on;
plot(y, x, 'r');
hold on;
plot(centerline_c, centerline_r, 'b.');
end