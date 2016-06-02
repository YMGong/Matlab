clear all;
close all;

fname=('/betaG/ControlOuter000000.2d.hdf5');
nlevel=2;

tarnumber=2;
ncomp=24;
tarcomp=tarnumber+1;

for level = 1:nlevel
    leveldata{level} = readchombolevel(fname,ncomp,level-1,tarcomp);
    [data{level}, data_mask{level}] = extractchombodata(fname,level,tarnumber); 
    nbox{level} = leveldata{level}.nbox;
    cd_i{level} = leveldata{level}(1).ii(1,2)-leveldata{level}(1).ii(1,1);
    cd_j{level} = leveldata{level}(1).jj(1,2)-leveldata{level}(1).jj(1,1);
    for ibox=1:nbox{level}
        y{level}{ibox}=leveldata{level}(ibox).ii(1:end);
        x{level}{ibox}=leveldata{level}(ibox).jj(1:end);
        numPbox_r{level}{ibox}=leveldata{level}(ibox).ni;
        numPbox_c{level}{ibox}=leveldata{level}(ibox).nj;
    end
    
end
if 0
for level= 1:nlevel
    data{level} = extractchombodata(fname,level,tarnumber); 
end

for level= 1:nlevel
    nbox{level} = leveldata{level}.nbox;
end

for level = 1:nlevel
   cd_i{level} = leveldata{level}(1).ii(1,2)-leveldata{level}(1).ii(1,1);
end

for level = 1:nlevel
   cd_j{level} = leveldata{level}(1).jj(1,2)-leveldata{level}(1).jj(1,1);
end
end

figure(1);
imagesc(data{1});hold on;
if 0
figure(2);
imagesc(data{2});hold on;
figure(3);
imagesc(data{3});hold on;
figure(4);
imagesc(data_mask{1});colorbar;
end

%if 0
figure(5);
for level=1:nlevel;
%level=2;
    for ibox=1:nbox{level}
        r_n=0;c_n=0;
        
        for i=1:length(data_mask{level}(:,1))
            for j=1:length(data_mask{level}(1,:))
                if data_mask{level}(i,j)==ibox
                    
                 r_n=i;c_n=j;
                 
                 
                else
                end
               
            end
        end
        d=data{level}(r_n-numPbox_r{level}{ibox}+1:r_n,c_n-numPbox_c{level}{ibox}+1:c_n);
        imagesc(x{level}{ibox},y{level}{ibox},d); hold on;
    
    end
    
    
end
axis equal tight;
%end
