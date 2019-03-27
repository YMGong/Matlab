function Basins=GenMask_Basins_SW_150m(Basins)
%Sidralik Br: 33910
Basins(Basins==Basins(16430, 3396) |...
    Basins==Basins(16480, 3573))=Basins(16570, 3356);

%Frederikshab: 47673

%Nakaissorssuaq: 20762
Basins(Basins==Basins(15420, 3037) |...
    Basins==Basins(15620, 3277))=Basins(15510, 2893);

%Kangiata nunata se:27660
Basins(Basins==Basins(14770, 4085) |...
       Basins==Basins(15090, 3500) )=Basins(14810, 2930);
temp=Basins(14530:14560, 2804:2898);
temp(temp==Basins(14550, 2839)) = Basins(14540, 2879);
Basins(14530:14560, 2804:2898)=temp;

%Narssap se:32993
Basins(Basins==Basins(14220, 3147) |...
    Basins==Basins(14420, 2982))=Basins(14200, 2925);

temp=Basins(14300:14410, 2879:3076);
temp(temp==Basins(14370, 2952)) = Basins(14200, 2925);
Basins(14300:14410, 2879:3076)=temp;

temp=Basins(14410:14430, 2923:2965);
temp(temp==Basins(14370, 2952)) = Basins(14200, 2925);
Basins(14410:14430, 2923:2965)=temp;

temp=Basins(14408:14430, 2927:2954);
temp(temp==Basins(14414, 3941)) = Basins(14200, 2925);
Basins(14408:14430, 2927:2954)=temp;

% %clean up
% temp=Basins(2492:2498, 461:471);
% temp(temp==1697) = 2203;
% Basins(2492:2498, 461:471)=temp;
end