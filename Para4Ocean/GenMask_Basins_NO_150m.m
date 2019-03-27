function Basins=GenMask_Basins_NO_150m(Basins)
%Humboldt Gl.:56675
% Basins(Basins==Basins(3349, 1761))=Basins(4233, 2893);
% Basins(Basins==Basins(3040, 1705))=Basins(4233, 2893);
% Basins(Basins==Basins(3401, 1960))=Basins(4233, 2893);
%Basins(Basins==Basins(2985, 1745))=Basins(4233, 2893);

Basins(Basins==Basins(3444, 2109))=Basins(4233, 2893);
Basins(Basins==Basins(2960, 1851))=Basins(4233, 2893);
Basins(Basins==Basins(2970, 1925))=Basins(4233, 2893);
Basins(Basins==Basins(2828, 1910))=Basins(4233, 2893); % the strap on the right
Basins(Basins==Basins(2818, 1976))=Basins(4233, 2893);

temp=Basins(2995:3070, 1672:1721);
temp(temp==Basins(4233, 2893)) = Basins(3057, 1691);
Basins(2995:3070, 1672:1721)=temp;

temp=Basins(2574:2735, 2189:2332);
temp(temp==Basins(4233, 2893)) = 0;
Basins(2574:2735, 2189:2332)=temp;

temp=Basins(2594:2787, 2334:2742);
temp(temp==Basins(4233, 2893)) = 0;
Basins(2594:2787, 2334:2742)=temp;
%Basins(Basins==Basins(2779, 1941))=Basins(4233, 2893);
%Basins(Basins==Basins(2735, 1894))=Basins(4233, 2893);

%Petermann Gl. :69480
temp=Basins(2523:2628, 2466:2608);
temp(temp==Basins(2558, 2562)) = Basins(2534, 2493);
Basins(2523:2628, 2466:2608)=temp;

Basins(Basins==Basins(2949, 2598))=Basins(3230, 2878);
Basins(Basins==Basins(2347, 2751))=Basins(3230, 2878);
Basins(Basins==Basins(2803, 3059))=Basins(3230, 2878);
temp=Basins(2236:2633, 2353:2857);
temp(temp==Basins(2415, 2589)) = Basins(3230, 2878);
Basins(2236:2633, 2353:2857)=temp;

temp=Basins(2255:2374, 2608:2744);
temp(temp==Basins(3230, 2878)) = Basins(2276, 2618);
Basins(2255:2374, 2608:2744)=temp;

temp=Basins(2460:2467, 2531:2537);
temp(temp==Basins(3230, 2878)) = Basins(2462, 2534);
Basins(2460:2467, 2531:2537)=temp;

%Steenby Gl.: 66017
Basins(Basins==Basins(1901, 3380))=Basins(2059, 3367);
Basins(Basins==Basins(2111, 3337))=Basins(2059, 3367);
Basins(Basins==Basins(2142, 3363))=Basins(2059, 3367);
temp=Basins(1853:1941, 3315:3386);
temp(temp~=Basins(3230, 2878) & temp~=0.0) = Basins(2059, 3367);
Basins(1853:1941, 3315:3386)=temp;

%Ryder Gl.: 17590
Basins(Basins==Basins(1880, 3804))=Basins(2388, 3908);

% C.H. Ostenfeld Gl.:57912
Basins(Basins==Basins(2195, 4279))=Basins(2300, 4673);
Basins(Basins==Basins(1917, 4481))=Basins(2300, 4673);
temp=Basins(1839:2218, 4212:4378);
temp(temp~=Basins(2300, 4673) & temp~=0.0) = Basins(2300, 4673);
Basins(1839:2218, 4212:4378)=temp;

% Hagen Gl. & Academy Gl. :18714
Basins(Basins==Basins(2043, 5451) |...
    Basins==Basins(2010, 5571) |...
    Basins==Basins(2068, 5848) |...
    Basins==Basins(2176, 5971) |...
    Basins==Basins(2010, 6091) |...
    Basins==Basins(2082, 6036) |...
    Basins==Basins(2253, 5610) )=Basins(1728, 6185);
%Basins==Basins(1737, 5667) |...
temp=Basins(1539:1633, 5532:5617);
temp(temp==Basins(1728, 6185)) = Basins(1588, 5533);
Basins(1539:1633, 5532:5617)=temp;

temp=Basins(2016:2075, 6166:6238);
temp(temp==Basins(1728, 6185)) = Basins(2020, 6230);
Basins(2016:2075, 6166:6238)=temp;

temp=Basins(1543:1820, 5886:5989);
temp(temp==Basins(1728, 6185)) = Basins(1652, 5908);
Basins(1543:1820, 5886:5989)=temp;

temp=Basins(1541:1783, 5942:6057);
temp(temp==Basins(1728, 6185)) = Basins(1652, 5908);
Basins(1541:1783, 5942:6057)=temp;

%Heilprin Gl. :18971
Basins(Basins==Basins(4361, 1329) |...
    Basins==Basins(4280, 1242))=Basins(4238, 1120);
temp=Basins(4186:4312, 1092:1259);
temp(temp~=Basins(4238, 1120) & temp~=0.0) = Basins(4238, 1120);
Basins(4186:4312, 1092:1259)=temp;

%Tracy Gl.:59220
temp=Basins(4114:4169, 1098:1169);
temp(temp~=Basins(4088, 1411) & temp~=0.0) = Basins(4088, 1411);
Basins(4114:4169, 1098:1169)=temp;
temp=Basins(4138:4226, 1158:1298);
temp(temp~=Basins(4088, 1411) & temp~=0.0) = Basins(4088, 1411);
Basins(4138:4226, 1158:1298)=temp;

%Farquhar Gl.:18972
Basins(Basins==Basins(4057, 1138))=Basins(4086, 1175);

temp=Basins(4080:4093, 1109:1122);
temp(temp==Basins(4089, 1115)) = Basins(4086, 1175);
Basins(4080:4093, 1109:1122)=temp;

temp=Basins(3810:4146, 1078:1433);
temp(temp==Basins(4086, 1175)) = 18972;
Basins(3810:4146, 1078:1433)=temp;

end