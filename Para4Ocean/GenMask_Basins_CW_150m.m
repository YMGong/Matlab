function Basins=GenMask_Basins_CW_150m(Basins)
% Jakobshavn Isbrae: 26614
Basins(Basins==Basins(10940, 3178) |...
    Basins==Basins(10770, 3399) |...
    Basins==Basins(10630, 3403))=Basins(11101, 3423);

% Sermeq K: 55774
Basins(Basins==Basins(10320, 3042) |...
    Basins==Basins(11340, 3507))=Basins(10300, 3250);

temp=Basins(10410:10560, 2950:3356);
temp(temp==Basins(10460, 3152)) = Basins(10170, 3079);
Basins(10410:10560, 2950:3356)=temp;

% Sermeq avangnardieq: 45493
Basins(Basins==Basins(10190, 3179))=Basins(10170, 3079);

temp=Basins(10245:10275, 2979:3049);
temp(temp~=Basins(10170, 3079) & temp~=0.0) = Basins(10170, 3079);
Basins(10245:10275, 2979:3049)=temp;

% Store: 9788
temp=Basins(9879:10050, 2927:3141);
temp(temp~=Basins(10170, 3079) & temp~=0.0) = Basins(9940, 3107);
Basins(9879:10050, 2927:3141)=temp;

temp=Basins(9871:10030, 2932:3048);
temp(temp==Basins(9971, 3002)) = Basins(9940, 3107);
Basins(9871:10030, 2932:3048)=temp;

%Rinks Isbrae: 56467
Basins(Basins==Basins(9101, 2943))=Basins(8755, 3450);

end