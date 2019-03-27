function Basins=mergeBasins(Basins, DB, mini_area)
flag1=1;
flag2=1;
nrDB = numel(unique(DB.Z(:)))-1; % nr of drainage basins
STATS = regionprops(DB.Z,'PixelIdxList','Area','Centroid');
for run = 1:nrDB
 if STATS(run).Area*DB.cellsize^2 < mini_area
     % get the color code and the coordinates of the small basins
    [x,y] = ind2coord(DB,...
        sub2ind(DB.size,...
        round(STATS(run).Centroid(2)),...
        round(STATS(run).Centroid(1))));
    small_patch(flag1,:)=[run;x;y];
    flag1=flag1+1;
 else
     % get the color code and the coordinates of the big basins
    [x,y] = ind2coord(DB,...
        sub2ind(DB.size,...
        round(STATS(run).Centroid(2)),...
        round(STATS(run).Centroid(1))));
    big_patch(flag2,:)=[run;x;y];
    flag2=flag2+1;
  end
end
% merge the small basins with the nearest big basins
for i = 1:length(small_patch)
    dist = sqrt((big_patch(:,2)-small_patch(i,2)).^2+...
        (big_patch(:,3)-small_patch(i,3)).^2);
    [~,index] = sort(dist);
    Basins(Basins==small_patch(i,1))=big_patch(index(1),1);
end
end