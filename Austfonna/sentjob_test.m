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
for stepV = 4:2:10; % stepping through the regular array for nearest
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

fprintf('saving files...\n')
filename=(['bed2d_trimmed_',num2str(figureID)]);
save(filename,'bed2d_trimmed');

filename=(['surf2d_',num2str(figureID+1)]);
save(filename,'surf2d');

filename=(['thick2d_trimmed_',num2str(figureID+2)]);
save(filename,'thick2d_trimmed');


figureID=figureID+3;
end
