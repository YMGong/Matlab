
close all;
clear all;

%********************friction & viscosity coefficient*********************%

fprintf('Reading data from hdf5 files...\n')

pathG='/wrk/ygong/BISICLES/ASF_BISICLES/';
pathS='/wrk/ygong/BISICLES/ASF_BISICLES/betaS/';
fname='ControlOuter000032.2d.hdf5';

ncomp=24;
level=0;
tarcomp=22; %the ordinal of the componant not the showing number!!!
Beta_G=readchombolevel([pathG,fname],ncomp,level,tarcomp);
Beta_S=readchombolevel([pathS,fname],ncomp,level,tarcomp);

fprintf('Displying the structure...\n')

Beta_G
Beta_S

fprintf('Plotting...\n')

figure(1);
subplot(2,1,1)
data=Beta_G.data;
imagesc(data);colorbar;set(gca,'YDir','normal');
subplot(2,1,2)
data=Beta_S.data;
imagesc(data);colorbar;set(gca,'YDir','normal');


