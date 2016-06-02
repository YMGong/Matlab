if 0
close all;
clear all;

fprintf('loading data...\n')
bed  = load('bed.xyz');
surf = load('surf.xyz');
vel_11_x = load('v95s_p11s.x');
vel_11_y = load('v95s_p11s.y');

% start and end locations of 33X in UTM coordinates
xs = 545000;  xe = 755000;
ys = 8800000; ye = 8950000;
stepD = 400; % stepping through the regular array for nearest
             % neighbour calculations (in practice this gives grid
             % size in metres, you probably want somewhere between
             % 100m - 500m for Austfonna)
stepV = 1; % stepping through the scattered data (set this to 1 to
           % use all the data, or a higher integer to skip some
           % data, which will speed up this script at the cost of
           % accuracy)

fprintf('interpolating bed onto regular 2d array...\n')
[bed2d, Xq, Yq]  = interp2array(xs,xe,ys,ye,bed,stepV,stepD,-2000.);

fprintf('interpolating surf onto regular 2d array...\n')
[surf2d, Xq, Yq] = interp2array(xs,xe,ys,ye,surf,stepV,stepD,-2000.);

fprintf('interpolating velocity onto regular 2d array...\n')
[vel_x_2d, Xq, Yq]  = interp2array(xs,xe,ys,ye,vel_11_x,stepV,stepD,-2000.);

[vel_y_2d, Xq, Yq] = interp2array(xs,xe,ys,ye,vel_11_y,stepV,stepD,-2000.);

Mvel = sqrt(vel_y_2d.^2+vel_x_2d.^2);
imagesc(Mvel);colorbar;

fprintf('trimming non-physical data for topography...\n')
mask2d = ones(size(bed2d));mask2d(bed2d>surf2d)=0.;
bed2d_trimmed = bed2d;bed2d_trimmed(mask2d==0.)=surf2d(mask2d==0.);
thick2d = surf2d - bed2d; 
thick2d_trimmed=thick2d;thick2d_trimmed(thick2d<0.)=NaN;

figure(1),
imagesc(surf2d);colorbar;


figure(2),
imagesc(vel_y_2d);colorbar;

end
fprintf('saving the original data...\n')

filename=('ASFdata_2011Ori_for_BISICLES.mat');
save(filename,'bed2d','surf2d','bed2d_trimmed','thick2d','thick2d_trimmed','vel_x_2d','vel_y_2d','Xq','Yq');

%figure(1);clf
%imagesc(Xq,Yq,bed2d); set(gca,'YDir','normal');colorbar
%mesh(Xq,Yq,bed2d);set(gca,'YDir','normal');colorbar
%hold on;
%plot3(bed(:,1),bed(:,2),bed(:,3),'o');

% saveas(gcf,['/wrk/ygong/BISICLES/',num2str(i),'.fig','fig']); 


%**************************sorting the data out***************************%
if 0
fprintf('sorting the data out...\n')

fprintf('generating mask...\n')

mini_thick=30;
thickness_mini_above=thick2d_trimmed;
thickness_mini_above(thick2d_trimmed<=mini_thick)=NaN;

mask_mirri = zeros(size(thickness_mini_above));
mask_mirri(~(isnan(thickness_mini_above))) = 1.;
mask_flip = flipud(mask_mirri);
mask = bwareaopen(mask_flip,5*10e3); % set the pixel for a approparite number;
                                   % Remove all objects in the image containing fewer than that pixels.
                                   % Remove VSF


mask(236:252,226:244) = 1;
mask(286:296,211:220) = 1;
mask(262:265,182:186) = 1;
mask(307:316,278:293) = 0;
mask(295:314,176:183) = 0;
mask(56:67,361:370) = 0;
mask(332:335,197:201) = 0;

fprintf('preparing data...\n')

thickness = flipud(thick2d_trimmed);
thickness(mask==0.) = 0.;

Vel_x=flipud(vel_x_2d);
Vel_x(mask==0.) = 0.;

Vel_y=flipud(vel_y_2d);
Vel_y(mask==0.) = 0.;

bed_topo=flipud(bed2d_trimmed);
bed_topo(mask==0.) = 0.;


%----------------------------------------------------------------------------------------------
fprintf('expanding data...\n')

row_exp=384; colume_exp=640; %a domain with dimensions that can be divided by 2 many times

ASFmask=mask;ASFmask(row_exp,colume_exp)=0;
topg=bed_topo;topg(row_exp,colume_exp)=0;
uvel=Vel_x;uvel(row_exp,colume_exp)=0;
vvel=Vel_y;vvel(row_exp,colume_exp)=0;
thk=thickness;thk(row_exp,colume_exp)=0;
lsrf=topg;
usrf=topg+thk;
figure(1);imagesc(ASFmask);colorbar;
filename=('ASFmask.mat');
save(filename,'ASFmask');
%beta_g = rho * g * H * |grad(s)| / |u| u is the velocity vector |u|the value
%of the vector
%beta_g = (1 + rhoi* grav * thck * sqrt(dsx**2 + dsy**2)) / ( max(1.0d-6,umod - umodsia))
%beta_g =100 at openocean
%umod = sqrt(uvel**2 + vvel**2) sqrt: kaifang
% flwa = 4.0e-17;umodsia = (2.0d0*flwa(1:ewn,1:nsn,5) *thck**(glen_n+1)) / (glen_n+1.0d0) * (rhoi* grav)**glen_n* sqrt(dsx**2 + dsy**2)**glen_n
%minbeta_g < beta_g < maxbeta_g  minbeta_g=20 maxbeta_g = 1.0e+4

rhoi=918.0;grav=9.81;
[dsx dsy]=gradient(usrf,200);
vlu_grdS=sqrt(dsx.^2+dsy.^2);
vlu_u=sqrt(uvel.^2+vvel.^2);
%beta_g=(1+rhoi*grav*thk.*vlu_grdS)./vlu_u;
beta_g=rhoi*grav*thk.*vlu_grdS./(max(1.0e-6,vlu_u));
%beta_g=rhoi*grav*thk.*vlu_grdS./vlu_u;

beta_g(thk<=0)=100;beta_g(beta_g<20)=20;beta_g(beta_g>1.0e+4)=1.0e+4;
%beta_g(thk<=0)=100;beta_g(beta_g<20)=20;beta_g((beta_g>1.0e+4)|(isnan(beta_g)))=1.0e+4;

usq = uvel.^2 + vvel.^2;

velc(1:row_exp,1:colume_exp)=0;
velc((thk > 0) & (usq > 10)) = 1;
%velc(row_exp,colume_exp)=0;

divuh(row_exp,colume_exp)=0;
divuhc(row_exp,colume_exp)=0;
%**************************************************************************
    

    beta_s(1:row_exp,1:colume_exp) = 3.0e+2;
    beta_s(((topg + thk) > 400.0)&(usq < 100.0))=1.0e+4;
   
%**************************************************************************

figure(2);imagesc(velc);colorbar
figure(3);imagesc(beta_g);colorbar
figure(4);imagesc(beta_s);colorbar


fprintf('saving data in ascii and mat format...\n')

save('lsrf_asf.dat','lsrf','-ascii' );
save('uvel_asf.dat','uvel','-ascii' );
save('vvel_asf.dat','vvel','-ascii' );
save('topg_asf.dat','topg','-ascii' );
save('usrf_asf.dat','usrf','-ascii' ); 


filename=('ASFfor_BISICLES.mat');
save(filename,'ASFmask','thk','uvel','vvel','topg');


%if 0
fprintf('creating netcdf files...\n')
fprintf('getting the data...\n')   
x=Xq;sigma=10;
for row= 1:length(Yq)
    y(row,1)=Yq(1,row);
end


nccreate('ASFfor_BISICLES.nc','x','Dimensions',{'x' colume_exp},'Format','classic');
ncwrite('ASFfor_BISICLES.nc','x',x);

nccreate('ASFfor_BISICLES.nc','y','Dimensions',{'y' row_exp},'Format','classic');
ncwrite('ASFfor_BISICLES.nc','y',y);

nccreate('ASFfor_BISICLES.nc','sigma','Dimensions',{'sigma' 10},'Format','classic');
ncwrite('ASFfor_BISICLES.nc','sigma',sigma);


%filename={bed_topo,'surf_topo','Vel_x','Vel_y','thickness_without_NaN','Vel_x_without_NaN','Vel_y_without_NaN'};
%nccreate('ASFdata_for_BISICLES.nc','thickness','Dimensions',{'y' 751 'x' 1051},'Format','classic');
%ncwrite('ASFdata_for_BISICLES.nc','thickness',thickness);



for i=1:9 % the number of the additional filenames
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
        case 3
            data=uvel;filename='xvel'; 
        case 4
            data=vvel;filename='yvel';
        case 5
            data=velc;filename='velc';   
        case 6
            data=divuh;filename='divuh';       
        case 7
            data=divuhc;filename='divuhc';  
        case 8
            data=beta_g;filename='beta_g';   
        case 9
            data=beta_s;filename='beta_s';  
            
    end
    nccreate('ASFfor_BISICLES.nc',filename,'Dimensions',{'y' row_exp 'x' colume_exp},'Format','classic');
    ncwrite('ASFfor_BISICLES.nc',filename,data);
    
end


ncdisp('ASFfor_BISICLES.nc');
end
fprintf('end\n')

















