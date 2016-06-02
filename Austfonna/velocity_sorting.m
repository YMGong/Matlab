%close all;
%clear all;

fprintf('loading data...\n');
load('velocity_interpolate.mat');
load('thick2d_trimmed_3.mat');

fprintf('setting nans in vel to be 0...\n');

%figure(1);
%vel_x_nan2o(:,:)=vel_x_2d;
%vel_x_nan2o(isnan(vel_x_2d)) = 0.;
%imagesc(vel_x_nan2o);colorbar;

%figure(2);
%vel_y_nan2o(:,:)=vel_y_2d;
%vel_y_nan2o(isnan(vel_y_2d)) = 0.;
%imagesc(vel_y_nan2o);colorbar;

fprintf('matching nan2o data with thickness...\n');

figure(3);
vel_x_nan2o_thickn(:,:)=vel_x_nan2o;
vel_x_nan2o_thickn(isnan(thickness_cleaned)) = 0.;
imagesc(vel_x_nan2o_thickn);colorbar;

figure(4);
vel_y_nan2o_thickn(:,:)=vel_y_nan2o;
vel_y_nan2o_thickn(isnan(thickness_cleaned)) = 0.;
imagesc(vel_y_nan2o_thickn);colorbar;

fprintf('matching original vel data with thickness nans...\n');

figure(5);
vel_x_thickn(:,:)=vel_x_2d;
vel_x_thickn(isnan(thickness_cleaned)) = nan;
imagesc(vel_x_thickn);colorbar;

figure(6);
vel_y_thickn(:,:)=vel_y_2d;
vel_y_thickn(isnan(thickness_cleaned)) = nan;
imagesc(vel_y_thickn);colorbar;

fprintf('matching original vel data with thickness nans and minimum thickness...\n');

minithickness=10;

figure(7);
vel_x_mini(:,:)=vel_x_2d;
vel_x_mini(isnan(thickness_cleaned)|(thickness_cleaned<=minithickness)) = nan;
imagesc(vel_x_mini);colorbar;

figure(8);
vel_y_mini(:,:)=vel_y_2d;
vel_y_mini(isnan(thickness_cleaned)|(thickness_cleaned<=minithickness)) = nan;
imagesc(vel_y_mini);colorbar;

fprintf('saving data...\n');
filename=('velocity_Sorted.mat');
save(filename,'vel_x_nan2o','vel_y_nan2o','vel_x_nan2o_thickn','vel_y_nan2o_thickn','vel_x_thickn','vel_y_thickn','vel_x_mini','vel_y_mini');

fprintf('end \n');
 