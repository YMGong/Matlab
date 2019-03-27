function savenc(outfilename, lon, lat, data, dataname)
    row = size(data,1);
    colume = size(data,2);
    nccreate([outfilename,'.nc'],'lon','Dimensions',{'x' row 'y' colume},'Format','netcdf4');
    
    ncwrite([outfilename,'.nc'],'lon',lon);
    ncwriteatt([outfilename,'.nc'],'lon','units','degree_east');
    ncwriteatt([outfilename,'.nc'],'lon','_CoordinateAxisType', "Lon");

    nccreate([outfilename,'.nc'],'lat','Dimensions',{'x' row 'y' colume},'Format','netcdf4');
    ncwrite([outfilename,'.nc'],'lat',lat);
    ncwriteatt([outfilename,'.nc'],'lat','units','degree_north');
    ncwriteatt([outfilename,'.nc'],'lat','_CoordinateAxisType', "Lat");
    
    nccreate([outfilename,'.nc'],dataname,'Dimensions',{'x' row 'y' colume},'Format','netcdf4');
    ncwrite([outfilename,'.nc'],dataname, data);
    ncwriteatt([outfilename,'.nc'],dataname,'coordinates', 'lat lon');


end