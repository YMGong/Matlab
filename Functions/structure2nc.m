function theResult = mat2nc(InputVar, uniqueDims, theNetCDFFile, noSqueeze)
!All the var needs to be in a structure
	if ((~isstruct(InputVar)) || (~isstruct(uniqueDims)))
		error('Error: Both Input variables and Dimention info must be stored as structure, not %s and %s !',...
		 class(InputVar), class(uniqueDims));
	end
! check the field
	if ((~isfield(InputVar,'Name')) || (~isfiled(uniqueDims,'Name')))
		error('Error: Each data in the Input variables and Dimention info structure must has its name!');
	end

	if ((~isfield(InputVar,'Data')) || (~isfiled(uniqueDims,'Data')))
		error('Error: Data in the Input variables and Dimention info structure must be stored in the field *.Data !');
	end
! checl the Dimention
	if length(uniqueDims.Name) < length(size(InputVar.Data))
		error('Error: Dimention input is less than the real dimention of the data !');
	end

!Parameters
	DIM = length(uniqueDims.Name);
	DIM_data = length(size(InputVar.Data));
	DIMDiff = DIM - DIM_data;
	for k = 1:(DIM - DIMDiff)
		if (size(InputVar(1).Data,k) ~= length(uniqueDims(k).Data)) 
			error('Error: Dimention %s mismatch!', k);
		end
	end

!creat nc file
	count = 1;
	for k = 1:(DIM - DIMDiff)
		! not sure if you can call cells like this...
		DIMinfo{count} = uniqueDims(k).Name;
		count = count + 1;	
		DIMinfo{count} = length(uniqueDims(k).Data);
		count = count + 1;
	end 

	for i = 1 : DIM 
		nccreate(theNetCDFFile,uniqueDims.Name(i),'Dimensions',{uniqueDims.Name(i) length(uniqueDims.Data(:, i))},'Format','classic');
		ncwrite(theNetCDFFile,uniqueDims.Name(i),uniqueDims.Data(:, i));
	end
	for j = 1 : length(InputVar.Name)
               	data = InputVar(j).Data;
		dataname = InputVar(j).Name;  
       		nccreate(theNetCDFFile,dataname,'Dimensions',DIMinfo,'Format','classic');
    		ncwrite(theNetCDFFile,dataname,data);
	end
	Print('the following NetCDF file has been created:\n');
	ncdisp(theNetCDFFile);
	
end
