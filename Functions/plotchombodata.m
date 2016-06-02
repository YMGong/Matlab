

function [] = plotchombodata(fname,data,nlevel,ncomp,tarnumber)

fprintf('data needs to be a structure contains the data in each level seperatly\n')
index=input('if data used is the same one contained in the file input 1; otherwise input 2 :\n');

tarcomp=tarnumber+1;
for level = 1:nlevel
    leveldata{level} = readchombolevel(fname,ncomp,level-1,tarcomp);
    [dataori{level}, data_mask{level}] = extractchombodata(fname,level,tarnumber); 
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



switch index
	   case 1
        data= dataori;  
       case 2
end

for level = 1:nlevel
figure(level);
imagesc(data{level});hold on;colorbar;title(['level',num2str(level)]);
end



    figure(level+1);
for level=1:nlevel;

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
colorbar;title('all levels');
fprintf('end\n');
end
