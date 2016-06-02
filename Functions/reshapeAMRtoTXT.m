
%*************************************************************************%
%**** Data needs to be created by the amrtotxt tool or at least in the ***%
%**** same format!!                                                    ***%
%*************************************************************************%
function [data_matrix] = reshapeAMRtoTXT(fname,var_index)

fprintf('data loading ... \n');
data = load(fname);

fprintf('data reshaping ... \n');

rows = unique(data(:,1),'rows');
i = size(rows,1);

columns = unique(data(:,2),'rows');
j = size(columns,1);

data_matrix = reshape(data(:,var_index),i,j);


end