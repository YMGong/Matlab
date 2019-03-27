function Basins=GenMask_Basins_NW_1km(Basins)
%Umiamak Isbrae :2299
%Ingia Isbrae: 803
temp=Basins(1406:1418, 458:467);
temp(temp==2299) = 803;
Basins(1406:1418, 458:467)=temp;
%Upernavik Isstrom: 265
Basins(Basins==228|Basins==412|Basins==372)=265;
temp=Basins(1358:1390, 434:461);
temp(temp==265) = 834;
Basins(1358:1390, 434:461)=temp;
%Unatakavsaup se: 458
Basins(Basins==212)=458;
%Kakivfait se: 857
temp=Basins(1223:1233, 392:411);
temp(temp==458) = 857;
Basins(1223:1233, 392:411)=temp;
temp=Basins(1215:1224, 464:504);
temp(temp==569) = 857;
Basins(1215:1224, 464:504)=temp;
%Qeqertarssup se: 569
temp=Basins(1178:1222, 355:502);
temp(temp==458|temp==249) = 569;
Basins(1178:1222, 355:502)=temp;
%249
%644
Basins(Basins==299)=644;
temp=Basins(1139:1147, 410:419);
temp(temp==249) = 644;
Basins(1139:1147, 410:419)=temp;

temp=Basins(1138:1143, 393:397);
temp(temp==249) = 644;
Basins(1138:1143, 393:397)=temp;

%1049
Basins(Basins==441|Basins==172|Basins==1049|Basins==213)=1049;
temp=Basins(1050:1056, 379:388);
temp(temp==645) = 1049;
Basins(1050:1056, 379:388)=temp;

%645
%466

%233
Basins(Basins==482|Basins==285)=233;
%444
%550
Basins(Basins==651)=550;
temp=Basins(890:911, 297:321);
temp(temp==366) = 550;
Basins(890:911, 297:321)=temp;
Basins(890, 322)=550;
%366
Basins(Basins==95|Basins==336)=366;
temp=Basins(882:898, 261:284);
temp(temp==406) = 366;
Basins(882:898, 261:284)=temp;
temp=Basins(888:898, 296:299);
temp(temp==550) = 366;
Basins(888:898, 296:299)=temp;

%406
%183
temp=Basins(871:877, 232:244);
temp(temp==470) = 183;
Basins(871:877, 232:244)=temp;

%470
Basins(Basins==266|Basins==330)=470;
end