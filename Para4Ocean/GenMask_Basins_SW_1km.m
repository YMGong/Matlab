function Basins=GenMask_Basins_SW_1km(Basins)
%Sidralik Br: 1432
Basins(Basins==1502)=1432;
%Uukkaasorsuaq: 1384
Basins(Basins==1514)=1384;
%Kangiata nunata se:284
Basins(Basins==685|Basins==1148)=284;
%Narssap se:1503
Basins(Basins==759|Basins==1860)=1503;

%clean up
temp=Basins(2492:2498, 461:471);
temp(temp==1697) = 2203;
Basins(2492:2498, 461:471)=temp;
end