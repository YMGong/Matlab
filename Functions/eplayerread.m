function [arrOut, Xq, Yq,data] = eplayerread(filename,l_s,l_e,colume,l_ss,l_se,s_colume,layers,data_index,mask)
% load ASF_B3_mask.mat;
% load ASF_synSMB.mat;
fprintf('Extrscting all the data from:\n');
fprintf([filename,'\n']);
fprintf('Coordinate...\n');
coord = readepdata(filename,l_s,l_e,colume,0);
fprintf('data...\n');
epdata = readepdata(filename,l_ss,l_se,s_colume,0);

[uniq_coord, ai, ci]=unique(coord,'rows');
data(:,1:3) = coord(:,1:3);
data(:,4) = epdata(:,data_index);
data_uniq(:,1:2) = data(ai,1:2);
data_uniq(:,3) = data(ai,4);


xs = 545000;  xe = 755000;
ys = 8800000; ye = 8950000;
stepD = 400;
stepV = 1;

data_layer(layers).data=0.0;
end
fprintf('extracting data and interpolating...\n');

for l = 1:layers
    fprintf(['layer: ',num2str(l),'\n']);
     data_raw=data_uniq(l:layers:end,:);
 
      [data_layer(l).data, Xq, Yq]  = interp2array(xs,xe,ys,ye,data_raw,stepV,stepD,-2000.,'natural');
     
      data_layer(l).data(row_exp,colume_exp) = 0.0; 
      data_layer(l).data(mask == 0)=0.0;
        
end

filename_mat = ([filename,'.mat']);
fprintf([filename,' has been saved.\n']);
save(filename_mat,'data_layer');
arrOut = data_layer;

imagesc(data_layer(11).data),colorbar;
 fprintf('end.\n')
end