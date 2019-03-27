function Basins=GenMask_Basins_NE_150m(Basins)
% 79North & Zachariae Isstrom:42495
Basins(Basins==Basins(2936, 6945) |...
    Basins==Basins(3013, 7312) |...
    Basins==Basins(3981, 7098) |...    
    Basins==Basins(3211, 7056))=Basins(2947, 7178);

temp=Basins(3288:3467, 7360:7440);
temp(temp==Basins(3394, 7401)) = Basins(2947, 7178);
Basins(3288:3467, 7360:7440)=temp;

temp=Basins(3344:3438, 7424:7471);
temp(temp==Basins(3394, 7401)) = Basins(2947, 7178);
Basins(3344:3438, 7424:7471)=temp;
% % Zachariae Isstrom :6397
% Basins(Basins==Basins(3544, 7336) |...
%     Basins==Basins(3409, 7409))=Basins(3692, 7310);
% 
% temp=Basins(2959:3463, 7345:7764);
% temp(temp~=Basins(3692, 7310) & temp~=0.0) = Basins(3692, 7310);
% Basins(2959:3463, 7345:7764)=temp;


% Storstrommen :9415
Basins(Basins==Basins(4903, 7277) |...
    Basins==Basins(4608, 7397) |...
    Basins==Basins(4424, 7520) |...
    Basins==Basins(4298, 7485) |...
    Basins==Basins(4157, 7656) |...
    Basins==Basins(4495, 7949) |...
    Basins==Basins(4143, 7885) |...
    Basins==Basins(4576, 6868) |...
    Basins==Basins(4790, 7114) )=Basins(4157, 7274);

temp=Basins(4557:4827, 7720:7923);
temp(temp==Basins(4157, 7274)) = Basins(4581, 7738);
Basins(4557:4827, 7720:7923)=temp;

temp=Basins(4465:4557, 7857:7895);
temp(temp==Basins(4157, 7274)) = Basins(4581, 7738);
Basins(4465:4557, 7857:7895)=temp;

temp=Basins(4579:4658, 7911:7940);
temp(temp==Basins(4157, 7274)) = Basins(4581, 7738);
Basins(4579:4658, 7911:7940)=temp;

temp=Basins(4517:4585, 7870:7919);
temp(temp==Basins(4157, 7274)) = Basins(4581, 7738);
Basins(4517:4585, 7870:7919)=temp;

temp=Basins(4255:4334, 7849:7993);
temp(temp==Basins(4581, 7738)) = Basins(4157, 7274);
Basins(4255:4334, 7849:7993)=temp;

temp=Basins(4711:4744, 7998:8078);
temp(temp==Basins(4157, 7274)) = Basins(4741, 8038);
Basins(4711:4744, 7998:8078)=temp;

temp=Basins(4245:4504, 7626:7916);
temp(temp==Basins(4504, 7916)) = Basins(4157, 7274);
Basins(4245:4504, 7626:7916)=temp;

temp=Basins(4495:4668, 7844:7986);
temp(temp==Basins(4504, 7916)) = Basins(4157, 7274);
Basins(4495:4668, 7844:7986)=temp;

temp=Basins(4477:4722, 7637:7971);
temp(temp==Basins(4558, 7754)) = Basins(4157, 7274);
Basins(4477:4722, 7637:7971)=temp;

% L.Bistrup Br. :18996
Basins(Basins==Basins(5309, 7986) |...
    Basins==Basins(5522, 7753) |...
    Basins==Basins(5224, 7841) |...
    Basins==Basins(5496, 7621) |...
    Basins==Basins(5602, 7300) |...
    Basins==Basins(5288, 7593) |...
    Basins==Basins(5192, 7485) |...
    Basins==Basins(5410, 7363) )=Basins(4952, 8106);

temp=Basins(5153:5252, 8166:8234);
temp(temp==Basins(4952, 8106)) = Basins(5254, 8178);
Basins(5153:5252, 8166:8234)=temp;

temp=Basins(5097:5220, 8210:8326);
temp(temp==Basins(4952, 8106)) = Basins(5254, 8178);
Basins(5097:5220, 8210:8326)=temp;

temp=Basins(5071:5123, 8228:8329);
temp(temp==Basins(4952, 8106)) = Basins(5254, 8178);
Basins(5071:5123, 8228:8329)=temp;

temp=Basins(4979:5119, 8265:8348);
temp(temp==Basins(4952, 8106)) = Basins(5254, 8178);
Basins(4979:5119, 8265:8348)=temp;

temp=Basins(5068:5151, 8227:8311);
temp(temp==Basins(4952, 8106)) = Basins(5254, 8178);
Basins(5068:5151, 8227:8311)=temp;

temp=Basins(5125:5159, 8200:8215);
temp(temp==Basins(4952, 8106)) = Basins(5159, 8215);
Basins(5125:5159, 8200:8215)=temp;

%Waltershausen Gl. :45311
Basins(Basins==Basins(6659, 7941) |...
    Basins==Basins(7938, 6362) |...
    Basins==Basins(6423, 7717) |...
    Basins==Basins(6298, 7900) |...
    Basins==Basins(6636, 7709) |...
    Basins==Basins(6806, 7800) |...
    Basins==Basins(6645, 7363) )=Basins(6640, 8149);

temp=Basins(6429:6447, 8084:8108);
temp(temp==Basins(6640, 8149)) = Basins(6433, 8101);
Basins(6429:6447, 8084:8108)=temp;

temp=Basins(6439:6464, 8184:8208);
temp(temp==Basins(6640, 8149)) = Basins(6433, 8101);
Basins(6439:6464, 8184:8208)=temp;

temp = Basins(6043:7321, 6685:8672);
temp(int64(temp)==45310) = 45311;
Basins(6043:7321, 6685:8672)=temp;

temp = Basins(6671:6757, 8255:8319);
temp(int64(temp)==45311) = 0;
Basins(6671:6757, 8255:8319)=temp;

temp = Basins(6680:6723, 8241:8257);
temp(int64(temp)==45311) = 0;
Basins(6680:6723, 8241:8257)=temp;

temp = Basins(6688:6714, 8220:8241);
temp(int64(temp)==45311) = 0;
Basins(6688:6714, 8220:8241)=temp;
end