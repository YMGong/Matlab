
%*************************************************************************%
% Using AMRtoTXT.sh to get the data from hdf5 files first                 %
%  Deleting the last line of the txt file by hand (need to find a way automaticelly)%
%*************************************************************************%

close all;
clear all;
%*************************************************************************%

fprintf('Loading data...\n');
filename = 'BISICLES_200_temp_t4_C';
%filename= 'BISICLES_400_temp_t4_C';
%filename1= 'BISICLES_400_Beta_t4_C';
filename1 = 'BISICLES_200_Beta_t4_C';
index_C = 3;
index_mu = 4;
%*************************************************************************%
% struName = [filename,'.name'];
% struTemp = [filename,'.temp'];
% struBeta = [filename1,'.beta'];
% struVelb = [filename1,'.velb'];
% struVelo = [filename1,'.velo'];
fprintf('Reshaping...\n');

layers = 12;
for l = 1:layers
    
    temp(l).name = ['layer = ',num2str(l),'%0.6g'];
    
    temp(l).temp = reshapeAMRtoTXT([filename,'.txt'],l+2);
    temp(l).temp = rot90(temp(l).temp);
    %tem_layer_old(l).temph = reshapeAMRtoTXT(filename,l+3);
    
end


figure(1),imagesc(temp(1).temp),colorbar;
figure(2),imagesc(temp(12).temp),colorbar;
figure(3),imagesc(temp(5).temp),colorbar;
eval([filename '= temp;']);% give the value of  temp to a variable has the name of filename

ctrl.beta = reshapeAMRtoTXT([filename1,'.txt'],index_C);
ctrl.beta = rot90(ctrl.beta);

ctrl.mu = reshapeAMRtoTXT([filename1,'.txt'],index_mu);
ctrl.mu = rot90(ctrl.mu);

eval([filename1 '= ctrl;']);
imagesc(ctrl.beta);colorbar;

fprintf('Saving...\n');
file_name = [filename,'.mat'];
save(file_name,'BISICLES_200_temp_t4_C');
 file_name1 = [filename1,'.mat'];
 save(file_name1,'BISICLES_200_Beta_t4_C');

fprintf('end.\n')
