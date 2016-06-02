delete Elmer2011velofor_BISICLES.nc
close all;
clear all;
load('ASFdata_2011Ori_for_BISICLES.mat');
load maskASF_Elmer.mat;
load maskB3.mat;
load ASFElmerdata.mat;
mask = Mask_ASF_Elmer;
%-----------------------------
fprintf('Elmer/Ice data loading...\n');
path1 = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/';
load ([path1,'asf_adjoint_ti_2011_6.ep.mat']);
load ([path1,'ASF_Elmer_400_t6_t5.mat']);
Elmer_2011Beta_t6_velu_surface = matrixdata(2).data_layer(11).data;
Elmer_2011Beta_t6_velu_surface(mask==0)=0.0;
Elmer_2011Beta_t6_velv_surface = matrixdata(3).data_layer(11).data;
Elmer_2011Beta_t6_velv_surface(mask==0)=0.0;
imagesc(sqrt(Elmer_2011Beta_t6_velu_surface.^2 + Elmer_2011Beta_t6_velv_surface.^2));
colorbar;colormap(jet);
%Elmer_2011Beta_t6_velw_surface = matrixdata(4).data_layer(11).data;
%-----------------------------
% mini_thick=30;
% thickness_mini_above=thick2d_trimmed;
% 
% thickness_mini_above(thick2d_trimmed<=mini_thick)=NaN;
% 
% mask_mirri = zeros(size(thickness_mini_above));
% mask_mirri(~(isnan(thickness_mini_above))) = 1.;
% mask_flip = flipud(mask_mirri);
% mask = bwareaopen(mask_flip,10e4); % set the pixel for a approparite number;
%                                       %Remove all objects in the image containing fewer than that pixels.
% 
% 
% mask(470:510,450:490) = 1;
% mask(570:590,420:440) =1;
% mask(613:632,572:584) =0;
% mask(606:624,352:362) =0;
% mask(520:532,352:364) =0;
% mask(110:130,723:738) =0;
% figure(1);imagesc(mask);colorbar;

%filename=('ASFmask.dat');
%save(filename,'mask','-ascii');
%mask = rot90(ASFmask,-1);




% fprintf('Saving data...\n')
% filename=('ASF2011for_BISICLES.mat');
% save(filename,'mask','thickness','Vel_x','Vel_y','bed_topo');
%-----------------------------------------------------------------------------------------------
fprintf('expanding data...\n')
obvely = Elmer_2011Beta_t6_velu_surface;
obvelx = Elmer_2011Beta_t6_velv_surface;
%row_exp=640; colume_exp=384; %a domain with dimensions that can be divided by 2 many times
row_exp=384; colume_exp=640;

topg = topo;
%surf = flipud(surf2d);
%surf = surf2d;surf(row_exp,colume_exp)=0;surf(mask==0.) = 0.;surf(isnan(surf)) = 0.;
if(isnan(surf))
    fprintf('there is nans in surf');
end
%topg=flipud(bed2d_trimmed);
%topg=bed2d_trimmed;topg(row_exp,colume_exp)=0;topg(mask==0.) = 0.;surf(isnan(surf)) = 0.;
if(isnan(topg))
    fprintf('there is nans in topg');
end
%uvel=flipud(vel_x_2d);
%uvel=vel_x_2d;uvel(row_exp,colume_exp)=0;uvel(mask==0.) = 0.;uvel(isnan(uvel)) = 0.;
if(isnan(obvelx))
    fprintf('there is nans in uvel');
end
%vvel=flipud(vel_y_2d);
%vvel=vel_y_2d;vvel(row_exp,colume_exp)=0;vvel(mask==0.) = 0.;vvel(isnan(vvel)) = 0.;
if(isnan(obvely))
    fprintf('there is nans in vvel');
end
%thk=flipud(thick2d_trimmed);
%mini_thick=30;
%thk=thick2d_trimmed;thk(row_exp,colume_exp)=0;thk(isnan(thk)) = 0.;thk(thk<mini_thick) = mini_thick;thk(mask==0.) = 0.;
if(isnan(thk))
    fprintf('there is nans in thk');
end

yvel_neg = -obvely;

lsrf=topg;
usrf=topg+thk;
%beta_g = rho * g * H * |grad(s)| / |u| u is the velocity vector |u|the value
%of the vector
%beta_g = (1 + rhoi* grav * thck * sqrt(dsx**2 + dsy**2)) / ( max(1.0d-6,vvel_asf.datumod - umodsia))
%beta_g =100 at openocean
%umod = sqrt(uvel**2 + vvel**2) sqrt: kaifang
% flwa = 4.0e-17;umodsia = (2.0d0*flwa(1:ewn,1:nsn,5) *thck**(glen_n+1)) / (glen_n+1.0d0) * (rhoi* grav)**glen_n* sqrt(dsx**2 + dsy**2)**glen_n
%minbeta_g < beta_g < maxbeta_g  minbeta_g=20 maxbeta_g = 1.0e+4

rhoi=910.0;grav=9.81;
[dsx dsy]=gradient(usrf,400);

%  dx_g(row_exp,colume_exp) = 0;
%  dx_g(:,1) = (usrf(:,2) - usrf(:,1))/200;
%  dx_g(:,colume_exp) = (usrf(:,colume_exp) - usrf(:,colume_exp-1))/200;
%  dx_g(:,2:colume_exp-1) = (usrf(:,3:colume_exp) - usrf(:,1:colume_exp-2))/400;

vlu_grdS=sqrt(dsx.^2+dsy.^2);
vlu_u=sqrt(obvelx.^2+obvely.^2);
%vlu_u=sqrt(velu95.^2+velv95.^2);
%make the interio stiker
beta_g=rhoi*grav*thk.*vlu_grdS./(max(1.0e-6,vlu_u-10));
beta_g(thk<=0)=100;beta_g(beta_g<20)=20;beta_g(beta_g>1.0e+5)=1.0e+5;

beta_gori=rhoi*grav*thk.*vlu_grdS./(max(1.0e-6,vlu_u));
beta_gori(thk<=0)=100;beta_gori(beta_gori<20)=20;beta_gori(beta_gori>1.0e+4)=1.0e+4;
%beta_g=rhoi*grav*thk.*vlu_grdS./vlu_u;

%beta_g(thk<=0)=100;beta_g(beta_g<20)=20;beta_g((beta_g>1.0e+4)|(isnan(beta_g)))=1.0e+4;
if(isnan(beta_g))
    fprintf('there is nans in beta_g');
end
%usq = uvel.^2 + vvel.^2;
usq = obvelx.^2 + obvely.^2;
Mvel = sqrt(usq);
%figure(2),imagesc(Mvel);colorbar;
velc(1:row_exp,1:colume_exp)=0;
velc((thk > 0) & (usq > 10)) = 1;
%velc(row_exp,colume_exp)=0;

divuh(row_exp,colume_exp)=0;
divuhc(row_exp,colume_exp)=0;
%**************************************************************************
    

    beta_elmer(1:row_exp,1:colume_exp) = 10^(-2) * 1.0e+6;
    beta_elmer(thk<=0)=1000;
   
%**************************************************************************

%**************************************************************************
    

    beta_s(1:row_exp,1:colume_exp) = 3.0e+2;
    beta_s(((topg + thk) > 400.0)&(usq < 100.0))=1.0e+4;
   
%************************************************************************vvel_asf.dat**
figure(1);imagesc(beta_g);colorbar;
figure(2);imagesc(beta_gori);colorbar;
%figure(3);imagesc(beta_g);colorbar;

%**************************************************************************


fprintf('preparing temperature data...\n')
temp_ocean = 272.0;
thk_pl = thk/10;
temp000000 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*0.5*thk_pl;
temp000000(mask==0.) = temp_ocean;temp000000(temp000000>temp_ocean) = temp_ocean;temp000000(isnan(temp000000)) = temp_ocean;
temp000001 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*1.5*thk_pl;
temp000001(mask==0.) = temp_ocean;temp000001(temp000001>temp_ocean) = temp_ocean;temp000001(isnan(temp000001)) = temp_ocean;
temp000002 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*2.5*thk_pl;
temp000002(mask==0.) = temp_ocean;temp000002(temp000002>temp_ocean) = temp_ocean;temp000002(isnan(temp000002)) = temp_ocean;
temp000003 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*3.5*thk_pl;
temp000003(mask==0.) = temp_ocean;temp000003(temp000003>temp_ocean) = temp_ocean;temp000002(isnan(temp000002)) = temp_ocean;temp000003(isnan(temp000003)) = temp_ocean;
temp000004 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*4.5*thk_pl;
temp000004(mask==0.) = temp_ocean;temp000004(temp000004>temp_ocean) = temp_ocean;temp000004(isnan(temp000004)) = temp_ocean;
temp000005 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*5.5*thk_pl;
temp000005(mask==0.) = temp_ocean;temp000005(temp000005>temp_ocean) = temp_ocean;temp000005(isnan(temp000005)) = temp_ocean;
temp000006 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*6.5*thk_pl;
temp000006(mask==0.) = temp_ocean;temp000006(temp000006>temp_ocean) = temp_ocean;temp000006(isnan(temp000006)) = temp_ocean;
temp000007 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*7.5*thk_pl;
temp000007(mask==0.) = temp_ocean;temp000007(temp000007>temp_ocean) = temp_ocean;temp000007(isnan(temp000007)) = temp_ocean;
temp000008 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*8.5*thk_pl;
temp000008(mask==0.) = temp_ocean;temp000008(temp000008>temp_ocean) = temp_ocean;temp000008(isnan(temp000008)) = temp_ocean;
temp000009 = 273.15-7.684 - 0.004*surf+40.0/2.072/1000*9.5*thk_pl;
temp000009(mask==0.) = temp_ocean;temp000009(temp000009>temp_ocean) = temp_ocean;temp000009(isnan(temp000009)) = temp_ocean;
if(isnan(temp000000))
    fprintf('there is nans in temp000000');
end
% figure(3);imagesc(temp000000);colorbar;
% figure(4);imagesc(temp000006);colorbar;
% figure(5);imagesc(temp000009);colorbar;


%**************************************************************************
%geothermal heat flux

geoheatflux(1:row_exp,1:colume_exp) = 0.04*365*24*60*60;
geoheatflux(mask==0.) = 0.;
%lambda(1:11).data(1:row_exp,1:colume_exp)=0;
for i=1:11
    lambda(i).data(1:row_exp,1:colume_exp) = mask;
    lambda(i).data(mask ==1) = 10^(i-1);
end
% lambda1 = mask;
% lambda1(mask == 1) = 1;
% lambda2 = mask;
% lambda2(mask == 1) = 10^2;
% lambda3 = mask;
% lambda3(mask == 1) = 10^3;
% lambda4 = mask;
% lambda4(mask == 1) = 10^4;
% lambda5 = mask;
% lambda5(mask == 1) = 10^5;
% lambda6 = mask;
% lambda6(mask == 1) = 10^6;
% lambda7 = mask;
% lambda7(mask == 1) = 10^7;
% lambda8 = mask;
% lambda8(mask == 1) = 108;
% lambda9 = mask;
% lambda9(mask == 1) = 10^9;
% lambda10 = mask;
% lambda10(mask == 1) = 10^10;
% figure(1);imagesc(thk);colorbar
% figure(2);imagesc(beta_g);colorbar
% figure(3);imagesc(beta_s);colorbar

% figure(6);imagesc(geoheatflux);colorbar;
% fprintf('saving data in ascii format...\n')
% 
% save('lsrf_asf.dat','lsrf','-ascii' );
% save('uvel_asf.dat','uvel','-ascii' );
% save('vvel_asf.dat','vvel','-ascii' );
% save('topg_asf.dat','topg','-ascii' );
% save('usrf_asf.dat','usrf','-ascii' ); 

fprintf('creating netcdf files...\n')
fprintf('getting the data...\n')   
y=Yq;y(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);
%x=rot90(Yq,-1);x(end:row_exp) = Yq(end,end):400:(Yq(end,end)+(row_exp-size(Yq,2))*400);
sigma=10;

x=rot90(Xq,-1);x(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);
%y=Xq;y(end:colume_exp) = Xq(end,end):400:(Xq(end,end)+(colume_exp-size(Xq,2))*400);


nccreate('Elmer2011velofor_BISICLES.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('Elmer2011velofor_BISICLES.nc','x',x);

nccreate('Elmer2011velofor_BISICLES.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('Elmer2011velofor_BISICLES.nc','y',y);

nccreate('Elmer2011velofor_BISICLES.nc','sigma','Dimensions',{'sigma' 10},'Format','classic');
ncwrite('Elmer2011velofor_BISICLES.nc','sigma',sigma);


%filename={bed_topo,'surf_topo','Vel_x','Vel_y','thickness_without_NaN','Vel_x_without_NaN','Vel_y_without_NaN'};
%nccreate('ASFdata_for_BISICLES.nc','thickness','Dimensions',{'y' 751 'x' 1051},'Format','classic');
%ncwrite('ASFdata_for_BISICLES.nc','thickness',thickness);

flag = 1;

for i=1:23 % the number of the additional filenames
    switch i
        
        case 1
            data=topg;filename='topg';
        %case 2
         %   data=surf_topo;filename='surf_topo';
        %case 3
         %   data=Vel_x;filename='Vel_x';
        %case 4
         %   data=Vel_y;filename='Vel_y';
        case 2
            data=thk;filename='thk';
%         case 3
%             data=uvel;filename='xvel'; 
%         case 4
%             data=vvel;filename='yvel';
        case 3
            data=velc;filename='velc';   
        case 4
            data=divuh;filename='divuh';       
        case 5
            data=divuhc;filename='divuhc';  
        case 6
            data=beta_g;filename='beta_g';   
        case 7
            data=beta_s;filename='beta_s';  
        case 8
            data=beta_elmer;filename='beta_elmer';  
        case 9
            data=temp000000;filename='temp000000';  
        case 10
            data=temp000001;filename='temp000001';  
        case 11
            data=temp000002;filename='temp000002';  
        case 12
            data=temp000003;filename='temp000003';  
        case 13
            data=temp000004;filename='temp000004';  
        case 14
            data=temp000005;filename='temp000005';  
        case 15
            data=temp000006;filename='temp000006';  
        case 16
            data=temp000007;filename='temp000007';  
        case 17
            data=temp000008;filename='temp000008';  
        case 18
            data=temp000009;filename='temp000009'; 
        case 19
            data=yvel_neg;filename='yvel_neg'; 
        case 20
            data=geoheatflux;filename='geoheatflux';
           
        case 21
            data=obvelx;filename='obvelx';
        case 22
            data=obvely;filename='obvely';   
        case 23
            data=beta_gori;filename='beta_gori';
    end
    nccreate('Elmer2011velofor_BISICLES.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp },'Format','classic');
    ncwrite('Elmer2011velofor_BISICLES.nc',filename,rot90(data,-1));
    flag = flag +1;
end
 for j=1:10
     data = lambda(j).data;filename=['lambda',num2str(j)];
     nccreate('Elmer2011velofor_BISICLES.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp },'Format','classic');
     ncwrite('Elmer2011velofor_BISICLES.nc',filename,rot90(data,-1));
 end

ncdisp('Elmer2011velofor_BISICLES.nc');

fprintf('end\n')






