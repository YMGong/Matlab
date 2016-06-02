%***************************Forward run results***************************%

clear all;
%close all;

path='/wrk/ygong/BISICLES/ASF_simulation/ASF_BISICLES_Forward_400/';

%fname='BISICLES.asf.800.001.hdf5';
fname='ASFForward_BISICELS_400.hdf5';
%h5disp([path,fname]);

level=0; % Always to be 0 after flattening
ncomp = h5readatt([path,fname],'/','num_components'); 

%need to read x_utm and y_utm out!

thickness_BISICLES_400_Nm = 0;
BISICLES_xVel_BISICLES_400_Nm = 1;
BISICLES_yVel_BISICLES_400_Nm = 2;
BISICLES_zVel_BISICLES_400_Nm = 3;
BISICLES_xfVel_BISICLES_400_Nm = 12;
BISICLES_yfVel_BISICLES_400_Nm = 13;
BISICLES_zfVel_BISICLES_400_Nm = 14;
BISICLES_xbVel_BISICLES_400_Nm = 15;
BISICLES_ybVel_BISICLES_400_Nm = 16;
BISICLES_zbVel_BISICLES_400_Nm = 17;
BISICLES_Z_surface_BISICLES_400_Nm = 4;
BISICLES_Z_bottom_BISICLES_400_Nm = 5;
BISICLES_Z_base_BISICLES_400_Nm = 6;
BISICLES_basal_friction_BISICLES_400_Nm = 7;

thickness_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,thickness_BISICLES_400_Nm+1);
BISICLES_xVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_xVel_BISICLES_400_Nm+1);
BISICLES_yVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_yVel_BISICLES_400_Nm+1);
BISICLES_zVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_zVel_BISICLES_400_Nm+1);
BISICLES_xfVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_xfVel_BISICLES_400_Nm+1);
BISICLES_yfVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_yfVel_BISICLES_400_Nm+1);
BISICLES_zfVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_zfVel_BISICLES_400_Nm+1);
BISICLES_xbVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_xbVel_BISICLES_400_Nm+1);
BISICLES_ybVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_ybVel_BISICLES_400_Nm+1);
BISICLES_zbVel_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_zbVel_BISICLES_400_Nm+1);
BISICLES_Z_surface_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_Z_surface_BISICLES_400_Nm+1);
BISICLES_Z_bottom_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_Z_bottom_BISICLES_400_Nm+1);
BISICLES_Z_base_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_Z_base_BISICLES_400_Nm+1);
BISICLES_basal_friction_BISICLES_400_LevelData = readchombolevel([path fname],ncomp,level,BISICLES_basal_friction_BISICLES_400_Nm+1);

dx = 400;
thickness_ini = flipud(thickness_BISICLES_400_LevelData.data);
BISICLES_xVel_ini = flipud(BISICLES_xVel_BISICLES_400_LevelData.data);
BISICLES_yVel_ini = flipud(BISICLES_yVel_BISICLES_400_LevelData.data);
BISICLES_zVel_ini = flipud(BISICLES_zVel_BISICLES_400_LevelData.data);
BISICLES_xfVel_ini = flipud(BISICLES_xfVel_BISICLES_400_LevelData.data);
BISICLES_yfVel_ini = flipud(BISICLES_yfVel_BISICLES_400_LevelData.data);
BISICLES_zfVel_ini = flipud(BISICLES_zfVel_BISICLES_400_LevelData.data);
BISICLES_xbVel_ini = flipud(BISICLES_xbVel_BISICLES_400_LevelData.data);
BISICLES_ybVel_ini = flipud(BISICLES_ybVel_BISICLES_400_LevelData.data);
BISICLES_zbVel_ini = flipud(BISICLES_zbVel_BISICLES_400_LevelData.data);
BISICLES_Z_surface_ini = flipud(BISICLES_Z_surface_BISICLES_400_LevelData.data);
BISICLES_Z_bottom_ini = flipud(BISICLES_Z_bottom_BISICLES_400_LevelData.data);
BISICLES_Z_base_ini = flipud(BISICLES_Z_base_BISICLES_400_LevelData.data);
BISICLES_basal_friction_ini = flipud(BISICLES_basal_friction_BISICLES_400_LevelData.data);

figure(2),imagesc(thickness_ini);
if 0
colume_exp = size(thickness_ini,1);
row_exp = size(thickness_ini,2);
fprintf('creating netcdf files...\n')
fprintf('getting the data...\n')   


%filename={bed_topo,'surf_topo','Vel_x','Vel_y','thickness_BISICLES_400_without_NaN','Vel_x_without_NaN','Vel_y_without_NaN'};
%nccreate('ASFdata_for_BISICLES.nc','thickness_ini','Dimensions',{'y' 751 'x' 1051},'Format','classic');
%ncwrite('ASFdata_for_BISICLES.nc','thickness_ini',thickness_ini);



for i=1:15% the number of the additional filenames
    switch i
        
        case 1
            data=thickness_ini;filename='thickness_ini';
        %case 2
         %   data=surf_topo;filename='surf_topo';
        %case 3
         %   data=Vel_x;filename='Vel_x';
        %case 4
         %   data=Vel_y;filename='Vel_y';
        case 2
            data=BISICLES_basal_friction_ini;filename='BISICLES_basal_friction_ini';
        case 3
            data=BISICLES_xVel_ini;filename='BISICLES_xVel_ini'; 
        case 4
            data=BISICLES_yVel_ini;filename='BISICLES_yVel_ini';
        case 5
            data=BISICLES_zVel_ini;filename='BISICLES_zVel_ini';   
        case 6
             data=BISICLES_xfVel_ini;filename='BISICLES_xfVel_ini';      
        case 7
             data=BISICLES_yfVel_ini;filename='BISICLES_yfVel_ini';
        case 8
             data=BISICLES_zfVel_ini;filename='BISICLES_zfVel_ini';   
        case 9
             data=BISICLES_xbVel_ini;filename='BISICLES_xbVel_ini';  
        case 10
             data=BISICLES_ybVel_ini;filename='BISICLES_ybVel_ini';
        case 11
             data=BISICLES_zbVel_ini;filename='BISICLES_zbVel_ini';   
        case 12
            data=BISICLES_Z_surface_ini;filename='BISICLES_Z_surface_ini';  
        case 13
            data=BISICLES_Z_bottom_ini;filename='BISICLES_Z_bottom_ini';  
        case 14
            data=BISICLES_Z_base_ini;filename='BISICLES_Z_base_ini';  
        case 15
            data=dx;filename='dx';
    end
    nccreate('ASFforward_BISICLES.nc',filename,'Dimensions',{'x' colume_exp 'y' row_exp },'Format','classic');
    ncwrite('ASFforward_BISICLES.nc',filename,data);
    
end

ncdisp('ASFforward_BISICLES.nc');
end

fprintf('saving...\n')
filename = ('ASFForward_BISICLES_400.mat');
save(filename,'thickness_ini','BISICLES_basal_friction_ini','BISICLES_xVel_ini','BISICLES_yVel_ini','BISICLES_zVel_ini','BISICLES_xfVel_ini','BISICLES_yfVel_ini','BISICLES_zfVel_ini','BISICLES_xbVel_ini','BISICLES_ybVel_ini','BISICLES_zbVel_ini','BISICLES_Z_surface_ini','BISICLES_Z_bottom_ini','BISICLES_Z_base_ini');

fprintf('end.\n')

