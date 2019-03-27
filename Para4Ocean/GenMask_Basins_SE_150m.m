function Basins=GenMask_Basins_SE_150m(Basins)
%Fenris Gl.: 26665
Basins(Basins==Basins(12650, 6619))=Basins(12430, 6496);

%Helheim Gl.:4677
Basins(Basins==Basins(12800, 6353))=Basins(12770, 5969);
temp=Basins(12880:12990, 6353:6471);
temp(temp~=Basins(12770, 5969) & temp~=0.0) = Basins(12770, 5969);
Basins(12880:12990, 6353:6471)=temp;

temp=Basins(12910:13100, 6255:6364);
temp(temp~=Basins(12770, 5969) & temp~=0.0) = Basins(12770, 5969);
Basins(12910:13100, 6255:6364)=temp;

temp=Basins(13040:13130, 6215:6304);
temp(temp~=Basins(12770, 5969) & temp~=0.0) = Basins(12770, 5969);
Basins(13040:13130, 6215:6304)=temp;

%Ikertivaq: 23816
Basins(Basins==Basins(13520, 5825) |...
    Basins==Basins(13360, 5610) |...
    Basins==Basins(13430, 5617) |...
    Basins==Basins(13560, 5666) )=Basins(13480, 5890);
temp=Basins(13589:13670, 5876:6007);
temp(temp~=Basins(13480, 5890) & temp~=0.0) = Basins(13480, 5890);
Basins(13589:13670, 5876:6007)=temp;

temp=Basins(13580:13620, 5953:6007);
temp(temp~=Basins(13480, 5890) & temp~=0.0) = Basins(13480, 5890);
Basins(13580:13620, 5953:6007)=temp;

%Koge Bugt Ikeq: 39056
Basins(Basins==Basins(13550, 5403))=Basins(13820, 5465);

temp=Basins(13760:14110, 5576:5736);
temp(temp==Basins(13900, 5621) & temp==Basins(14060, 5591)) = Basins(13820, 5465);
Basins(13760:14110, 5576:5736)=temp;

temp=Basins(13660:14060, 5473:5613);
temp(temp==Basins(13900, 5621) & temp==Basins(14060, 5591)) = Basins(13820, 5465);
Basins(13660:14060, 5473:5613)=temp;

temp=Basins(13751:14070, 5507:5733);
temp(temp==Basins(13910, 5652)) = Basins(13820, 5465);
Basins(13751:14070, 5507:5733)=temp;

temp=Basins(13670:14040, 5457:5591);
temp(temp==Basins(13910, 5652)) = Basins(13820, 5465);
Basins(13670:14040, 5457:5591)=temp;

%Gyldenlove: 63077
Basins(Basins==Basins(14480, 5323)|...
       Basins==Basins(14650, 5389)|...
       Basins==Basins(14750, 5258)|...
       Basins==Basins(14800, 4876)|...
       Basins==Basins(14800, 5530))=Basins(14480, 5323);
temp=Basins(14560:14620, 5404:5514);
temp(temp~=Basins(14480, 5323) & temp~=0.0) = Basins(14480, 5323);
Basins(14560:14620, 5404:5514)=temp;

%Bernstorff Isfj: 41360
Basins(Basins==Basins(14890, 5328))=Basins(14960, 5099);
temp=Basins(14930:14980, 5451:5476);
temp(temp~=Basins(14960, 5099) & temp~=0.0) = Basins(14960, 5099);
Basins(14930:14980, 5451:5476)=temp;

temp=Basins(14850:14880, 5447:5463);
temp(temp==Basins(14860, 5455)) = Basins(14960, 5099);
Basins(14850:14880, 5447:5463)=temp;

%Skinfaxe & Rimfaxe & GUldfaxe: 29970
Basins(Basins==Basins(15200, 4926))=Basins(15380, 5102);

%Heimdal Gl.: 11855
temp=Basins(15610:15830, 4957:5216);
temp(temp==Basins(15560, 5040)) = Basins(15540, 5006);
Basins(15610:15830, 4957:5216)=temp;

temp=Basins(15470:15620, 5005:5129);
temp(temp==Basins(15560, 5040)) = Basins(15540, 5006);
Basins(15470:15620, 5005:5129)=temp;

%Tingmiarmiut Fj. :25493
temp=Basins(15680:15820, 4861:4955);
temp(temp~=Basins(15750, 4876) & temp~=0.0) = Basins(15750, 4876);
Basins(15680:15820, 4861:4955)=temp;

temp=Basins(15730:15820, 4918:5004);
temp(temp~=Basins(15750, 4876) & temp~=0.0) = Basins(15750, 4876);
Basins(15730:15820, 4918:5004)=temp;

temp=Basins(15830:15850, 4932:4957);
temp(temp~=Basins(15750, 4876) & temp~=0.0) = Basins(15750, 4876);
Basins(15830:15850, 4932:4957)=temp;

% for Heimdal Gl.: 11855
temp=Basins(15590:15830, 4959:5223);
temp(temp==Basins(15560, 5040)) = Basins(15540, 5006);
Basins(15590:15830, 4959:5223)=temp;

end