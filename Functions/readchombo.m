function [a] = readchombo(fname)

% plot.pigv5.1km.l1l2.1lev.000159.2d.hdf5

% read number of variables and number levels in file
ncomp = h5readatt(fname,'/','num_components');    
nlevel = h5readatt(fname,'/','num_levels');

% read data level by level
for level = 1:nlevel
    %u{level} = readchombolevel(fname,ncomp,level-1,2);
    %v{level} = readchombolevel(fname,ncomp,level-1,3);
    Leveldata_G{level} = readchombolevel(fname,ncomp,level-1,1);    
end

% plotting
for level = 1:nlevel
    
    nbox = Leveldata_G{level}(1).nbox;
    
    for ibox = 1:nbox
        
    x = Leveldata_G{level}(ibox).ii(2:end-1);
    y = Leveldata_G{level}(ibox).jj(2:end-1);
    d = Leveldata_G{level}(ibox).data(2:end-1,2:end-1);
    
    imagesc(x,y,d); hold on; 
    
    end
   if 0
    for ibox = 1:nbox
    
    x = a{level}(ibox).i(2:end-1,2:end-1);
    y = a{level}(ibox).j(2:end-1,2:end-1);
    uu = u{level}(ibox).data(2:end-1,2:end-1);
    vv = v{level}(ibox).data(2:end-1,2:end-1);

    quiver(x,y,uu,vv,0,'k'); hold on; 
    end
    end
 end
if 0 
for level = 1:nlevel
    
    nbox = a{level}(1).nbox; 
    
    c = pickcolor(double(level-1)/2+0.25);
    
    for ibox = 1:nbox
            x = a{level}(ibox).ii;
            y = a{level}(ibox).jj;
            xl = x(1) + a{level}(ibox).dx/2;
            xu = x(end) - a{level}(ibox).dx/2;
            yl = y(1) + a{level}(ibox).dx/2;
            yu = y(end)- a{level}(ibox).dx/2;


            plot([xl xl xu xu xl],[yl yu yu yl yl],'Color',c)
    end
   
end
end
axis equal tight
colorbar

hold off

     %surf(a(ibox).i,a(ibox).j,sq(a(ibox).data(1,:,:))); hold on; pause