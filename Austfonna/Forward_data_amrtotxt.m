%*************************************************************************%
% Using AMRtoTXT.sh to get the data from hdf5 files first                 %
%  Deleting the last line of the txt file by hand (need to find a way automaticelly)%
%*************************************************************************%

close all;
clear all;
%*************************************************************************%

fprintf('Loading data...\n');

Elmer_forw(1).name='ASFForward_Elmer_500.txt';
Elmer_forw(2).name='ASFForward_Elmer_1500.txt';
Elmer_forw(3).name='ASFForward_Elmer_500_pl.txt';
Elmer_forw(4).name='ASFForward_Elmer_1500_pl.txt';

%*************************************************************************%

fprintf('Reshaping...\n');

xbvel_index = 3;
ybvel_index = 4;
zbvel_index = 5;
thickness_index = 6;

for i=1:4
Elmer_forw(i).xvel = reshapeAMRtoTXT(Elmer_forw(i).name,xbvel_index);
Elmer_forw(i).yvel = reshapeAMRtoTXT(Elmer_forw(i).name,ybvel_index);
Elmer_forw(i).zvel = reshapeAMRtoTXT(Elmer_forw(i).name,zbvel_index);
Elmer_forw(i).thickness = reshapeAMRtoTXT(Elmer_forw(i).name,thickness_index);
end

figure(1),imagesc(Elmer_forw(1).xvel),colorbar;
figure(2),imagesc(Elmer_forw(2).zvel),colorbar;
figure(3),imagesc(Elmer_forw(4).thickness),colorbar;

fprintf('Saving...\n');
filename = 'ASFForward_Elmer.mat';
save(filename,'Elmer_forw');
fprintf('end.\n')