function Basins=GenMask_Basins_CW_1km(Basins)
% Alangordliup se: 988
temp=Basins(1769:1796, 499:526);
temp(temp==568) = 988;
Basins(1769:1796, 499:526)=temp;
% Jakobshavn Isbrae: 568
Basins(Basins==186|Basins==112)=568;
temp=Basins(1737:1747, 536:542);
temp(temp==1318) = 568;
Basins(1737:1747, 536:542)=temp;
% Sermeq avangnardieq: 1318
% Sermeq K & a: 801
Basins(Basins==215|Basins==439|Basins==215)=801;
% Store: 666
Basins(Basins==1258)=666;
temp=Basins(1600:1610, 511:516);
temp(temp==865|temp==439) = 666;
Basins(1600:1610, 511:516)=temp;
% Sermeq s : 865
%Kangerdluarssup se :440
temp=Basins(1488:1546, 484:544);
temp(temp==865) = 440;
Basins(1488:1546, 484:544)=temp;
%kangerdlugssup se: 1086
temp=Basins(1455:1501, 491:530);
temp(temp==2299) = 1086;
Basins(1455:1501, 491:530)=temp;
%Rinks Isbrae: 541
temp=Basins(1444:1459, 489:506);
temp(temp==2299) = 541;
Basins(1444:1459, 489:506)=temp;

end