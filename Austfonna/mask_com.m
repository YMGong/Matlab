if 0
close all;
clear all;

fprintf('loading data...\n')
mask_R = load('ASFiceclean3.txt');
load('ASFmask.mat');
bed  = load('bed.xyz');
surf = load('surf.xyz');

xs = 545000;  xe = 755000;
ys = 8800000; ye = 8950000;
stepD = 200; % stepping through the regular array for nearest
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
  
fprintf('trimming non-physical data for topography...\n')
mask2d = ones(size(bed2d));mask2d(bed2d>surf2d)=0.;
bed2d_trimmed = bed2d;bed2d_trimmed(mask2d==0.)=surf2d(mask2d==0.);
thick2d = surf2d - bed2d; 
thick2d_trimmed=thick2d;thick2d_trimmed(thick2d<0.)=NaN;
           
fprintf('interpolating bed onto regular 2d array...\n')
[Mask, Xq, Yq]  = interp2array(xs,xe,ys,ye,mask_R,stepV,stepD,-2000.);
end
figure(1);
subplot(1,2,1);
imagesc(flipud(thick2d_trimmed));title('thick2d_trimmed');hold on;
subplot(1,2,2);
plot(mask_R(:,1),mask_R(:,2),'o');title('mask_R');hold on;