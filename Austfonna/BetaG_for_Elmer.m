
close all;
clear all;
load('ASFdata_2011Ori_for_BISICLES.mat');


lsrf=bed2d_trimmed;
usrf=bed2d_trimmed+thick2d_trimmed;
rhoi=910.0;grav=9.81;
[dsx dsy]=gradient(usrf,200);

vlu_grdS=sqrt(dsx.^2+dsy.^2);
vlu_u=sqrt(vel_x_2d.^2+vel_y_2d.^2);
%beta_g=(1+rhoi*grav*thk.*vlu_grdS)./vlu_u;
beta_g(:,:,1)=rhoi*grav*thick2d_trimmed.*vlu_grdS./(max(1.0e-6,vlu_u));
%beta_g=rhoi*grav*thk.*vlu_grdS./vlu_u;

beta_g(thick2d_trimmed<=0)=100;beta_g(beta_g<20)=20;beta_g(beta_g>1.0e+4)=1.0e+4;
beta_ori = beta_g;
beta_g=rot90(beta_g,-1);
%imagesc(beta_g);colormap;

t = 1; time =1;
var_t = linspace(1,t,t);
Yq(1,1:end) = Yq(1,end:-1:1); 

[y_utm,x_utm] = meshgrid(Yq,Xq);
x_utm=double(x_utm);y_utm=double(y_utm);
x = size(beta_g,1);
y = size(beta_g,2);
% nccreate('ASF_initBeta_forElmer.nc','x','Dimensions',{'x' x},'Format','classic');
% ncwrite('ASF_initBeta_forElmer.nc','x',x);
% 
% nccreate('ASF_initBeta_forElmer.nc','y','Dimensions',{'y' y},'Format','classic');
% ncwrite('ASF_initBeta_forElmer.nc','y',y);
% 
% nccreate('ASF_initBeta_forElmer.nc','t','Dimensions',{'t' t},'Format','classic');
% ncwrite('ASF_initBeta_forElmer.nc','t',t);
% 
% 
% %filename={bed_topo,'surf_topo','Vel_x','Vel_y','thickness_without_NaN','Vel_x_without_NaN','Vel_y_without_NaN'};
% %nccreate('ASFdata_for_BISICLES.nc','thickness','Dimensions',{'y' 751 'x' 1051},'Format','classic');
% %ncwrite('ASFdata_for_BISICLES.nc','thickness',thickness);
% 
% 
% 
% for i=1:3 % the number of the additional filenames
%     switch i
%         
% %         case 1
% %             data=t;filename='t_n';
% %         case 2
% %             data=time;filename='time';
%         case 1
%             data=x_utm;filename='x_utm'; 
%         case 2
%             data=y_utm;filename='y_utm';
%         case 3
%             data=beta_g;filename='beta_init';
%             
%     end
%     nccreate('ASF_initBeta_forElmer.nc',filename,'Dimensions',{'y' y 'x' x},'Format','classic');
%     ncwrite('ASF_initBeta_forElmer.nc',filename,data);
%     
% end
% nccreate('ASF_initBeta_forElmer.nc','t_n','Dimensions',{'y' 1 'x' 1},'Format','classic');
% ncwrite('ASF_initBeta_forElmer.nc','t_n',t);
% 
% nccreate('ASF_initBeta_forElmer.nc','time','Dimensions',{'y' 1 'x' 1},'Format','classic');
% ncwrite('ASF_initBeta_forElmer.nc','time',time);
figure(1);imagesc(beta_g);colorbar;
figure(2);imagesc(x_utm);colorbar;title('x_utm');
figure(3);imagesc(y_utm);colorbar;title('y_utm');


ncName = ('ASF_initBeta_forElmer.nc');
ncid   = netcdf.create(ncName,'clobber');
dimidx = netcdf.defDim(ncid,'x',x) ;
dimidy = netcdf.defDim(ncid,'y',y) ;
dimidt = netcdf.defDim(ncid,'t',t) ;
%dimidt = netcdf.defDim(ncid,'t',netcdf.getConstant('UNLIMITED')) ;
varidx = netcdf.defVar(ncid,'x_utm',netcdf.getConstant('double'),[dimidx,dimidy]) ;
varidy = netcdf.defVar(ncid,'y_utm',netcdf.getConstant('double'),[dimidx,dimidy]) ;
varidt2= netcdf.defVar(ncid,'t',netcdf.getConstant('int'),[dimidt]) ;
varidt = netcdf.defVar(ncid,'time',netcdf.getConstant('int'),[dimidt]) ;
varid  = netcdf.defVar(ncid,'beta_ini',netcdf.getConstant('double'),[dimidx,dimidy,dimidt]) ;

netcdf.endDef(ncid);

netcdf.putVar(ncid,varidx,x_utm);
netcdf.putVar(ncid,varidy,y_utm);
netcdf.putVar(ncid,varidt,var_t);
netcdf.putVar(ncid,varidt2,var_t);
netcdf.putVar(ncid,varid,beta_g);
netcdf.reDef(ncid);
netcdf.close(ncid);

ncdisp('ASF_initBeta_forElmer.nc');



fprintf('end\n')
