function Basins=GenMask_Basins_CE_1km(Basins)
%Daugaard Gl. :336
Basins(Basins==2681)=336;
temp=Basins(1361:1375, 1270:1280);
temp(temp==467) = 336;
Basins(1361:1375, 1270:1280)=temp;
% Rolige Br. :1838
Basins(Basins==2661)=1838;
%Vestford :611
Basins(Basins==1096|Basins==2777)=611;
temp=Basins(1533:1538, 1304:1309);
temp(temp==1838) = 611;
Basins(1533:1538, 1304:1309)=temp;
%Magga Dan Gl. :1387
Basins(Basins==333|Basins==2824)=1387;

%King Chrictian and many: 2238
Basins(Basins==1765)=2238;
%Kangerlugssuaq :389
Basins(Basins==244)=389;
%Polaric Gl. :562
Basins(Basins==1707)=562;
%Hutchinson plateau:303
temp=Basins(1870:1893, 1142:1225);
temp(temp==277) = 303;
Basins(1870:1893, 1142:1225)=temp;

%Deception:409
temp=Basins(1894:1901, 1210:1218);
temp(temp==277) = 409;
Basins(1894:1901, 1210:1218)=temp;

temp=Basins(1894:1906, 1170:1183);
temp(temp==277) = 409;
Basins(1894:1906, 1170:1183)=temp;

temp=Basins(1875:1886, 1157:1187);
temp(temp==409) = 303;
Basins(1875:1886, 1157:1187)=temp;

temp=Basins(1874:1880, 1150:1158);
temp(temp==409) = 303;
Basins(1874:1880, 1150:1158)=temp;

%Kruuse Fj. :277
temp=Basins(1940:1956, 1130:1152);
temp(temp==777) = 277;
Basins(1940:1956, 1130:1152)=temp;

temp=Basins(1953:1966, 1145:1161);
temp(temp==777) = 277;
Basins(1953:1966, 1145:1161)=temp;

Basins(Basins==227)=277;

%K.J.V. Steenstrup & Midgard Gl. :777

temp=Basins(1941:1958, 1130:1132);
temp(temp==277) = 777;
Basins(1941:1958, 1130:1132)=temp;

end