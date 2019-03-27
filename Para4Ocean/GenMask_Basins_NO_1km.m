function Basins=GenMask_Basins_NO_1km(Basins)
%Humboldt Gl.:173
Basins(Basins==142)=173;
temp=Basins(545:551, 332:339);
temp(temp==173) = 132;
Basins(545:551, 332:339)=temp;

%Petermann Gl. :340
temp=Basins(460:476, 444:473);
temp(temp==253) = 340;
Basins(460:476, 444:473)=temp;
Basins(Basins==264)=340;
%Steenby Gl.: 830

%Ryder Gl.:433
Basins(Basins==558)=433;
temp=Basins(377:426, 622:654);
temp((temp~=433)&(temp>0)) = 433;
Basins(377:426, 622:654)=temp;

% Harder Gl.:159
Basins(Basins==72)=159;
temp=Basins(364:370, 710:726);
temp(temp==159) = 936;
Basins(364:370, 710:726)=temp;
% Marie Sophie Gl.:2045
Basins(Basins==2017)=2045;
% Academy Gl. :94
temp=Basins(350:393, 899:921);
temp(temp==1752) = 1360;
Basins(350:393, 899:921)=temp;
Basins(Basins==1360)=94;
% Hagen Gl.:1752
Basins(Basins==1422)=1752;
Basins(Basins==913)=1752;
Basins(Basins==1857)=1752;
temp=Basins(388:394, 892:900);
temp(temp==1752) = 2045;
Basins(388:394, 892:900)=temp;
%Farquhar Gl.:18971
Basins(Basins==1422)=Basins(388:394, 892:900);

%Tracy Gl.:398
temp=Basins(721:756, 234:262);
temp(temp==117) = 398;
Basins(721:756, 234:262)=temp;
%Heilprin Gl. :1712
temp=Basins(735:745, 230:245);
temp(temp==117) = 996;
Basins(735:745, 230:245)=temp;
Basins(Basins==996)=1712;
temp=Basins(734:755, 230:255);
temp(temp==398) = 1712;
Basins(734:755, 230:255)=temp;
temp=Basins(745:757, 253:261);
temp(temp==398) = 1712;
Basins(745:757, 253:261)=temp;
end