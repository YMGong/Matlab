function Basins=GenMask_Basins_CE_150m(Basins)
%Daugaard Gl. :45310
Basins(Basins==Basins(8526, 8025) |...
    Basins==Basins(8581, 7659) |...
    Basins==Basins(8521, 7608) |...
    Basins==Basins(8391, 7754) |...
    Basins==Basins(8407, 7822) |...
    Basins==Basins(9269, 6617))=Basins(8644, 7670);
temp=Basins(8354:8542, 7902:8176);
temp(temp==Basins(8417, 8067)) = Basins(8644, 7670);
Basins(8354:8542, 7902:8176)=temp;

temp=Basins(8351:8359, 7942:7967);
temp(temp==Basins(8417, 8067)) = Basins(8352, 7951);
Basins(8351:8359, 7942:7967)=temp;

temp=Basins(8214:8368, 7876:8017);
temp(temp==Basins(8260, 7995)) = Basins(8644, 7670);
Basins(8214:8368, 7876:8017)=temp;

temp=Basins(8889:8923, 7805:7838);
temp(temp==Basins(8898, 7822)) = Basins(8644, 7670);
Basins(8889:8923, 7805:7838)=temp;

%Vestford :70318
%use the faster outlet
Basins(Basins==Basins(9904, 8057) |...   
    Basins==Basins(9829, 8194) |... 
    Basins==Basins(9611, 8019) |... 
    Basins==Basins(9527, 7994) |... 
    Basins==Basins(9768, 8004))=Basins(9664, 7831);

temp=Basins(9473:9672, 7968:8307);
temp(temp~=Basins(9664, 7831) & temp~=0.0) = Basins(9664, 7831);
Basins(9473:9672, 7968:8307)=temp;

%King Chrictian: 44893
Basins(Basins==Basins(10841, 8211) |...   
    Basins==Basins(10811, 8305) |... 
    Basins==Basins(10561, 8279) |... 
    Basins==Basins(10671, 8034) |... 
    Basins==Basins(9768, 8004))=Basins(11011, 8351);
temp=Basins(11070:11169, 8290:8425);
temp(temp==Basins(11011, 8351)) = Basins(11011, 8351);
Basins(11070:11169, 8290:8425)=temp;

temp=Basins(10030:10410, 7512:8343);
temp(temp==Basins(10200, 7874)|temp==Basins(10140, 8103)) = Basins(11070, 8290);
Basins(10030:10410, 7512:8343)=temp;

%Kangerlugssuaq :45484
temp=Basins(11065:11090, 7585:7633);
temp(temp==Basins(11075, 7619)) = Basins(10980, 6917);
Basins(11065:11090, 7585:7633)=temp;

temp=Basins(11054:11070, 7562:7594);
temp(temp==Basins(11075, 7619)) = Basins(10980, 6917);
Basins(11054:11070, 7562:7594)=temp;


temp=Basins(11020:11140, 7254:7451);
temp(temp~=Basins(10980, 6917) & temp~=0.0) = Basins(10980, 6917);
Basins(11020:11140, 7254:7451)=temp;

%Deception:13995
Basins(Basins==Basins(11960, 7182) |...  
    Basins==Basins(11770, 7239) |... 
    Basins==Basins(11460, 7025))=Basins(12030, 7508);

temp=Basins(11660:11780, 7317:7385);
temp(temp~=Basins(12030, 7508) & temp~=0.0) = Basins(12030, 7508);
Basins(11660:11780, 7317:7385)=temp;

%Kruuse Fj. :20527
temp=Basins(12100:12210, 7261:7669);
temp(temp~=Basins(12150, 7092) & temp~=0.0) = Basins(12150, 7092);
Basins(12100:12210, 7261:7669)=temp;

%K.J.V. Steenstrup Nordre & Sonder Gl. :38464
Basins(Basins==Basins(12450, 7384) |...  
    Basins==Basins(12270, 6906))=Basins(12550, 7213);


end