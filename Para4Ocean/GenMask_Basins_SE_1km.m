function Basins=GenMask_Basins_SE_1km(Basins)
%Midgard Gl.:399

%Fenris Gl.: 1506
temp=Basins(2017:2047, 1045:1062);
temp(temp==399) = 1506;
Basins(2017:2047, 1045:1062)=temp;

temp=Basins(1985:2007, 1056:1072);
temp(temp==399) = 1506;
Basins(1985:2007, 1056:1072)=temp;

%Helheim Gl.:344
temp=Basins(2041:2059, 1025:1035);
temp(temp==972) = 344;
Basins(2041:2059, 1025:1035)=temp;
%Ikertivaq
Basins(Basins==670|Basins==602|Basins==522|...
    Basins==757|Basins==1273|Basins==1344)=796;
%Koge Bugt Ikeq: 1559
Basins(Basins==695|Basins==392|Basins==567)=1559;

temp=Basins(2149:2185, 891:908);
temp(temp==796) = 1559;
Basins(2149:2185, 891:908)=temp;
temp=Basins(2170:2194, 908:926);
temp(temp==796) = 1559;
Basins(2170:2194, 908:926)=temp;
%Gyldenlove:
Basins(Basins==304)=434;

%Bernstorff Isfj: 598

%Rimfaxe : 579
temp=Basins(2405:2410, 863:884);
temp(temp==574) = 579;
Basins(2405:2410, 863:884)=temp;

% GUldfaxe: 574

%Heimdal Gl.: 318
temp=Basins(2432:2462, 803:840);
temp(temp==574) = 318;
Basins(2432:2462, 803:840)=temp;

temp=Basins(2451:2463, 840:844);
temp(temp==574) = 318;
Basins(2451:2463, 840:844)=temp;


%Tingmiarmiut Fj. :380

%Mogens Heinesen Fj. :424
Basins(Basins==195|Basins==435)=424;
%282
%420
end