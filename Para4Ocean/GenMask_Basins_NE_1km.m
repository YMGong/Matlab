function Basins=GenMask_Basins_NE_1km(Basins)
% 79North :368
Basins(Basins==663)=368;
Basins(Basins==793|Basins==351)=368;
% Zachariae Isstrom :326
temp=Basins(546:650,1170:1226);
temp(temp==368) = 326;
Basins(546:650,1170:1226)=temp;
% Storstrommen :243
Basins(Basins==1613|Basins==1735)=243;
temp=Basins(815:826,1247:1255);
temp(temp==243) = 1316;
Basins(815:826,1247:1255)=temp;

temp=Basins(756:766,1194:1206);
temp(temp==243) = 243;
Basins(756:766,1194:1206)=temp;

temp=Basins(767:776,1202:1213);
temp(temp==243) = 243;
Basins(767:776,1202:1213)=temp;

% L.Bistrup Br. :2480
temp=Basins(811:855, 1259:1293);
temp(temp==243) = 2480;
Basins(811:855, 1259:1293)=temp;
temp=Basins(830:843, 1194:1202);
temp(temp==2480) = 2517;
Basins(830:843, 1194:1202)=temp;
Basins(Basins==2490|Basins==2971)=2480;
%Waltershausen Gl. :2914
Basins(Basins==2683|Basins==3104)=2914;

%clean up


end