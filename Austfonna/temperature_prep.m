close all;
clear all;
load('ASFdataOri_for_BISICLES.mat');
load ASFmask.mat;
thk = ncread('ASFfor_BISICLES.nc','thk');
mask = ASFmask;


fprintf('expanding data...\n')

row_exp=384; colume_exp=640; %a domain with dimensions that can be divided by 2 many times
surf_el = flipud(surf2d);
surf_el(row_exp,colume_exp)=0;
surf_el(mask==0.) = 0.;
figure(1);imagesc(surf_el);colorbar;
figure(2);imagesc(thk);colorbar;

fprintf('preparing temperature data...\n')

temp000000 = -7.684 - 0.004*surf_el;