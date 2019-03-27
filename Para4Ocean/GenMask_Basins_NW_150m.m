function Basins=GenMask_Basins_NW_150m(Basins)
%Upernavik Isstrom: 41515
Basins(Basins==Basins(8219, 2462) |...
    Basins==Basins(8165, 2632) |...
    Basins==Basins(8452, 2905) |...
    Basins==Basins(8486, 2555) |...
    Basins==Basins(7888, 2553))=Basins(8096, 3212);

temp=Basins(7982:8112, 2301:2440);
temp(temp~=Basins(8096, 3212) & temp~=0.0) = Basins(8096, 3212);
Basins(7982:8112, 2301:2440)=temp;

%Unatakavsaup se: 71925
Basins(Basins==Basins(7855, 2283))=Basins(7743, 2326);

%K: 361
Basins(Basins==Basins(7447, 2288) |...
    Basins==Basins(7337, 2296))=Basins(7743, 2326);

temp=Basins(7507:7697, 2119:2416);
temp(temp==Basins(7649, 2285)) = Basins(7743, 2326);
Basins(7507:7697, 2119:2416)=temp;

temp=Basins(7162:7447, 2118:2952);
temp(temp==Basins(7319, 2453)) = Basins(7743, 2326);
Basins(7162:7447, 2118:2952)=temp;

temp=Basins(7382:7646, 2092:2387);
temp(temp==Basins(7319, 2453)) = Basins(7743, 2326);
Basins(7382:7646, 2092:2387)=temp;

%: 50283
Basins(Basins==Basins(7168, 2394) |...
    Basins==Basins(7015, 2238))=Basins(6900, 3092);

%25089
Basins(Basins==Basins(6719, 2349) |...
    Basins==Basins(6564, 2209) |...
    Basins==Basins(6283, 2202) |...
    Basins==Basins(6348, 2038) |...
    Basins==Basins(6177, 2174) |...
    Basins==Basins(6037, 2318))=Basins(6054, 3309);
%17146
Basins(Basins==Basins(5993, 2129) |...
    Basins==Basins(5756, 1898) |...
    Basins==Basins(5698, 2187))=Basins(5709, 2638);

%6700
Basins(Basins==Basins(5452, 1898))=Basins(5400, 2126);
temp=Basins(5602:5701, 1768:1961);
temp(temp==Basins(5640, 1842)) = Basins(5400, 2126);
Basins(5602:5701, 1768:1961)=temp;
% 
% 60835
Basins(Basins==Basins(5507, 1734))=Basins(5392, 1793);

%33151
temp=Basins(5405:5438, 1615:1665);
temp(temp==Basins(5431, 1649) | temp==Basins(5413, 1634)) = Basins(5274, 1795);
Basins(5405:5438, 1615:1665)=temp;

temp=Basins(5408:5425, 1605:1625);
temp(temp==Basins(5413, 1634)) = Basins(5274, 1795);
Basins(5408:5425, 1605:1625)=temp;

%31427
Basins(Basins==Basins(5047, 1437))=Basins(5081, 1689);
temp=Basins(5088:5377, 1287:1538);
temp(temp==Basins(5182, 1457)) = Basins(5081, 1689);
Basins(5088:5377, 1287:1538)=temp;

temp=Basins(4996:5262, 1473:1623);
temp(temp==Basins(5182, 1457)) = Basins(5081, 1689);
Basins(4996:5262, 1473:1623)=temp;

temp=Basins(5011:5244, 1527:1625);
temp(temp==Basins(5147, 1569)) = Basins(5081, 1689);
Basins(5011:5244, 1527:1625)=temp;

%57168
Basins(Basins==Basins(5068, 1815) |...
    Basins==Basins(5327, 1610))=Basins(5095, 1875);
temp=Basins(5255:5426, 1541:1671);
temp(temp==Basins(5338, 1606)) = Basins(5095, 1875);
Basins(5255:5426, 1541:1671)=temp;

temp=Basins(5072:5343, 1607:1857);
temp(temp==Basins(5338, 1606)) = Basins(5095, 1875);
Basins(5072:5343, 1607:1857)=temp;
% 

% %72367
% temp=Basins(5074:5147, 1231:1278);
% temp(temp==Basins(5103, 1238)) = Basins(4943, 1273);
% Basins(5074:5147, 1231:1278)=temp;
% 
% %56299
% 
% %44553
% temp=Basins(5103:5208, 1023:1075);
% temp(temp==Basins(5140, 1049)) = Basins(4986, 1032);
% Basins(5103:5208, 1023:1075)=temp;
% 
% %3009
% Basins(Basins==Basins(4894, 927))=Basins(4989, 860);
% %54246
% %34248
% temp=Basins(4712:4863, 562:684);
% temp(temp==Basins(4792, 615)) = Basins(4885, 574);
% Basins(4712:4863, 562:684)=temp;
% 
% %39000

end