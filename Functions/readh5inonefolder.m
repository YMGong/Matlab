function data = readh5inonefolder(fdir,fname_front,fname_back,level,tarnumber,fnumb,order_n,num_str,num_end)

fprintf('If the filename is not with consequent number, the first file must be the NetCDF file and all the NetCDF files must be line up together!\n')
index=input('if the filename is continuous with consequent number input 1; otherwise input 2 :\n');

strtl=order_n + 2; endl=order_n + 2 + fnumb - 1;

switch index
	   case 1
           fprintf('case 1\n');  
	         flag=1;
		     for i=num_str:num_end
			    data(:,:,flag)=extractchombodata([fname_front,num2str(i),fname_back],level,tarnumber);
		        flag=flag+1;
             end

        case 2
         fprintf('case 2\n'),file_stru=dir(fdir);  
		
       		for j=strtl:endl
               		
	            filename=file_stru(j).name;
                f=extractchombodata(filename,level,tarnumber);
                flag=j-3;
                data(:,:,flag)=f(:,:);
		    end

end

end