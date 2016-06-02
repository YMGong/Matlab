function [data, data_mask] = extractchombodata(fname,level,tarnumber)

ncomp = h5readatt(fname,'/','num_components');    
tarcomp=tarnumber+1;
Leveldata_G = readchombolevel(fname,ncomp,level-1,tarcomp);   

  nbox = Leveldata_G(1).nbox;
 
 
  CD_r=Leveldata_G(1).ii(1,2)-Leveldata_G(1).ii(1,1);
  CD_c=Leveldata_G(1).jj(1,2)-Leveldata_G(1).jj(1,1);
  
  
  baselevel =  readchombolevel(fname,ncomp,0,tarcomp);  

  base_start_r=baselevel(1).i(1,1);
  base_start_c=baselevel(1).j(1,1);
  
  


for ibox = 1:nbox
    % Finding the starting position of row in the mask
    x_s = Leveldata_G(ibox).ii(1,1);
    posi_rs = (x_s-base_start_r)/CD_r+1;
    % Finding the ending position of row in the mask
    x_e = Leveldata_G(ibox).ii(1,end);
    posi_re = (x_e-base_start_r)/CD_r+1;
    
    % Finding the starting position of row in the mask
    y_s = Leveldata_G(ibox).jj(1,1);
    posi_cs = (y_s-base_start_c)/CD_c+1;
    % Finding the ending position of row in the mask
    y_e = Leveldata_G(ibox).jj(1,end);
    posi_ce = (y_e-base_start_c)/CD_c+1;
    baselevel(1).nbox;
    d = Leveldata_G(ibox).data(1:end,1:end);
    
    x = Leveldata_G(ibox).ii(1:end);% Plotting without line 
    y = Leveldata_G(ibox).jj(1:end);
    %figure(1);imagesc(x,y,d); hold on; colorbar
    
    data(posi_rs:posi_re,posi_cs:posi_ce)=flipud(rot90(d));
   %data(posi_rs:posi_re,posi_cs:posi_ce)=d;
    data_mask(posi_rs:posi_re,posi_cs:posi_ce)=ibox;
    %mask=flipud(rot90(data_mask));
end


%figure(1);
%imagesc(data);colorbar;


end