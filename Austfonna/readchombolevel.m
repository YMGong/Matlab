function a = readchombolevel(fname,ncomp,level,tarcomp)

% read all boxes' data from a given level 

% construct name of level
levelname = ['/level_' num2str(level) '/'];

% read grid spacing for this level
dx = h5readatt(fname,levelname,'dx');
 
% read information on boxes on this level
% this return their coordinates of upper and lower edges in x,y
% units are number of grid cells from origin
box = h5read(fname,[levelname 'boxes']);

% how many boxes are there
nbox = size(box.lo_i,1);

% loop through boxes
for ibox = 1:nbox

% disp([box.lo_i(ibox)-1 box.hi_i(ibox)+1 box.lo_j(ibox)-1 box.hi_j(ibox)+1])

% find sizes of box from upper and lower coords PLUS one ghost cell on each edge    
    a(ibox).ni = (box.hi_i(ibox) - box.lo_i(ibox) + 1);
    a(ibox).nj = (box.hi_j(ibox) - box.lo_j(ibox) + 1);

% find total number of points in box    
    a(ibox).size = a(ibox).ni * a(ibox).nj;
    
%     a(ibox).ii = [dx/2.0:dx:dx/2.0 + (a(ibox).ni-1) * dx];
%     a(ibox).jj = [dx/2.0:dx:dx/2.0 + (a(ibox).nj-1) * dx];
    
% create 1d array with coords of cell centres
    a(ibox).ii = dx * [box.lo_i(ibox):box.hi_i(ibox)] + dx/2;
    a(ibox).jj = dx * [box.lo_j(ibox):box.hi_j(ibox)] + dx/2;

% mesh grid may be useful but not at moment   
    [a(ibox).i,a(ibox).j] = meshgrid(a(ibox).ii,a(ibox).jj); 
    
end
 
% now read pointers to where data is held in 1d array (one pointer per box)
offs = int32(h5read(fname,[levelname 'data:offsets=0'])); 

% read the data itself in 1d array
data = h5read(fname,[levelname 'data:datatype=0']);     

en = 0;

% visit each box and extract data
for ibox = 1:nbox

% find pointer to start of box's data    
     en = offs(ibox);
     
% extract each variable in turn stepping forward using known number of points in box     
     for icomp = 1:ncomp
         
         st = en + 1; en = en + a(ibox).size;

% keep only requested variable         
% need to remember to change ordering of data ji -> ij         

         if icomp == tarcomp
            a(ibox).data(:,:) = reshape(data(st:en),a(ibox).ni,a(ibox).nj)';
         end
                 
     end

% also store this level's grid spacing and the number of boxes present     
     a(ibox).nbox = nbox;
     a(ibox).dx = dx;
 
end
 

 
 %hold on