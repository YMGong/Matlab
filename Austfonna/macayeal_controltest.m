close all;
clear all; 
% In Elmer/Ice the units are always MPa-yr-m
delete BISICLE_MacAyealTest.nc
ncfilename = 'BISICLE_MacAyealTest.nc';

yearinsec = 365.25*24*60*60; 
%for meshing
Lx = 200.0e3;
elementx = 50;
Ly = 50.0e3;
elementy = 25;
[Coordsx,Coordsy] = meshgrid(0:Lx/elementx:Lx,Ly:-Ly/elementy:0);
sigma = 2;
%for vertical coord
zs=500.0-1.0e-3*Coordsx+20.0*(sin(3.0*pi*Coordsx/Lx).*sin(2.0*pi*Coordsy/Ly));
zb=zs-1500.0+2.0e-3*Coordsx;
thk=zs-zb;
%for beta initial value
F1=sin(3.0*pi*Coordsx/Lx).*sin(pi*Coordsy/Ly);
F2=sin(pi*Coordsx/(2.0*Lx)).*cos(4.0*pi*Coordsy/Ly);
beta=5.0e3*F1+5.0e3*F2;
betaSquare=beta.*beta/yearinsec;
betaIni(1:size(betaSquare,1),1:size(betaSquare,2)) = 1.0e3/sqrt(yearinsec);

%for velocity
% velox = 4.753e-6*yearinsec*(sin(2.0*pi*(Ly-Coordsy)/Ly)+2.5*sin(pi*(Ly-Coordsy)/Ly));
% veloy = 0.0;
% figure,imagesc(betaIni),colorbar,colormap(jet);
% figure,imagesc(zs),colorbar,colormap(jet);
% figure,imagesc(zb),colorbar,colormap(jet);
% figure,imagesc(velo),colorbar,colormap(jet);

%data creating for bisicles
structuredata(1).name = 'Coordsx';structuredata(1).data = Coordsx;
structuredata(2).name = 'Coordsy';structuredata(2).data = Coordsy;
structuredata(3).name = 'zs';structuredata(3).data = zs;
structuredata(4).name = 'zb';structuredata(4).data = zb;
structuredata(5).name = 'thk';structuredata(5).data = thk;
structuredata(6).name = 'betaSquare';structuredata(6).data = betaSquare;
structuredata(7).name = 'betaIni';structuredata(7).data = betaIni;

nccreat4BISICLES(0:Lx/elementx:Lx,Ly:-Ly/elementy:0,sigma,ncfilename,structuredata);




