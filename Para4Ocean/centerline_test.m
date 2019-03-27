clear all;
close all;
inum = 20;
jnum = 20;
data = nan(inum, jnum);
n_thresh = 20;
step = 1;
%% assign some mask area
count = 0;

for j = 3:1:11
    i = 10;
    data(max(i-count,1):min(i+count,inum), j) = 1;   
    count = count+1;
end
count = 0;
for j = 19:-1:11
    i = 10;
    data(max(i-count,1):min(i+count,inum), j) = 1;   
    count = count+1;
end
%% the starting point and ending point
i_s = 5; i_e = 19; %y
j_s = 14; j_e = 8; %x
pt1 = [i_s, j_s];
pt2 = [i_e, j_e];

n = max(j_e -j_s,n_thresh);
t = linspace(0,1,n+1); 
%t = t(2:(end-1)); 
t = t(1:end); 
v = pt2 - pt1;
x = pt1(1) + t*v(1);    % p(t) = p1 + t*(p2-p1)
y = pt1(2) + t*v(2);
figure, imagesc(data), hold on;
plot(y, x, 'r'); hold on;
%% find the perpendicular lines
width = 1;

for r= 2:n
    v = width*v / norm(v);
    
    r_ll = max(round(x(r)+v(2)),1); 
    r_ll = min(r_ll, inum); 
    
    c_ll = max(round(y(r)-v(1)), 1);
    c_ll = min(c_ll, jnum);
    
    r_ur = max(round(x(r)-v(2)), 1);
    r_ur = min(r_ur, inum);
    
    c_ur = max(round(y(r)+v(1)), 1);
    c_ur = min(c_ur, jnum);
    
    while (~isnan( data( r_ll, c_ll ) ) || ...
             ~isnan( data( r_ur, c_ur ) ) )
         
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
        c_new = c_ll:step:c_ur;
        r_new = round(slope * (0:step:(c_ur - c_ll)) + r_ll); 
    else
        c_new = c_ll:-step:c_ur;
        r_new = round(slope * (0:-step:(c_ur - c_ll)) + r_ll); 
    end
    counter = 1;
    for k = 1:length(c_new)
        if ~isnan(data(r_new(k), round(c_new(k)))) 
            data(r_new(k), round(c_new(k))) = 2;
            r_tmp(counter) = r_new(k);
            c_tmp(counter) = round(c_new(k));
            counter = counter + 1;
        end
    end
    data(r_tmp(end) + round((r_tmp(1) - r_tmp(end)) / 2), c_tmp(1) + round((c_tmp(end) - c_tmp(1)) / 2)) = 3;
    clear r_tmp c_tmp
end
figure, h=imagesc(data); pbaspect([1 1 1]);colorbar;
set(h, 'alphadata', ~isnan(data));
hold on;
plot(y, x, 'r');