%  Sensetive test

close all;
clear all;

fprintf('loading data...\n')
bed  = load('bed.xyz');
surf = load('surf.xyz');

% start and end locations of 33X in UTM coordinates
xs = 545000;  xe = 755000;
ys = 8800000; ye = 8950000;
figureID=1;
for stepV = 4:2:20; % stepping through the regular array for nearest
             % neighbour calculations (in practice this gives grid
             % size in metres, you probably want somewhere between
             % 100m - 500m for Austfonna)
stepD = 200; % stepping through the scattered data (set this to 1 to
           % use all the data, or a higher integer to skip some
           % data, which will speed up this script at the cost of
           % accuracy)

fprintf('interpolating bed onto regular 2d array...\n')
[bed2d, Xq, Yq]  = interp2array(xs,xe,ys,ye,bed,stepV,stepD,-2000.);

fprintf('interpolating surf onto regular 2d array...\n')
[surf2d, Xq, Yq] = interp2array(xs,xe,ys,ye,surf,stepV,stepD,-2000.);

fprintf('trimming non-physical data...\n')
mask2d = ones(size(bed2d));mask2d(bed2d>surf2d)=0.;
bed2d_trimmed = bed2d;bed2d_trimmed(mask2d==0.)=surf2d(mask2d==0.);
thick2d = surf2d - bed2d; 
thick2d_trimmed=thick2d;thick2d_trimmed(thick2d<0.)=NaN;


figure(figureID);clf
mesh(Xq,Yq,bed2d);set(gca,'YDir','normal');colorbar
hold on;
plot3(bed(:,1),bed(:,2),bed(:,3),'o');
saveas(gcf,['/wrk/ygong/BISICLES/','stepV',num2str(figureID),'.fig'],'fig'); 

figure(figureID+1);clf
mesh(Xq,Yq,surf2d);set(gca,'YDir','normal');colorbar
hold on;
plot3(surf(:,1),surf(:,2),surf(:,3),'o');
saveas(gcf,['/wrk/ygong/BISICLES/','stepV',num2str(figureID+1),'.fig'],'fig'); 

figure(figureID+2);clf
subplot(2,2,1)
imagesc(Xq,Yq,bed2d); set(gca,'YDir','normal');colorbar
subplot(2,2,2)
imagesc(Xq,Yq,surf2d); set(gca,'YDir','normal');colorbar
subplot(2,2,3)
imagesc(Xq,Yq,thick2d); set(gca,'YDir','normal');colorbar
subplot(2,2,4)
imagesc(Xq,Yq,thick2d_trimmed); set(gca,'YDir','normal');colorbar
saveas(gcf,['/wrk/ygong/BISICLES/','stepV',num2str(figureID+2),'.fig'],'fig'); 


%for i=1:figureID+2
%     print(gcf,'-djpeg',['/wrk/ygong/BISICLES/',num2str(i),'.jpeg']); 
%end
% Saving image in fig--turn it on when needed
%for i=1:figureID+2
%     saveas(gcf,['/wrk/ygong/BISICLES/',num2str(i),'.fig'],'fig'); 
%end
figureID=figureID+3;
end
% Saving image in jgp--turn it on when needed













