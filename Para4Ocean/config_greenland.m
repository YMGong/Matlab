% constants for all the files
%% for polarstereo_fwd
EARTHRADIUS = 6378137.0;
ECCENTRICITY = 0.08181919;
LAT_TRUE = 70;
LON_POSY = -45;
%% ice properties
rhoi = 910.0;
rhow = 1028.0;
%% basin numbers
NO = [56675, 69480, 66017, 17590, 57912, 18714, 42495, 9415, 18996,];
% 18971, 59220, 18972 should belong to the north but the old dem in the
% region is collected from 1985
% 42495,  9415, 18996 should be in the NE but the old dem there is
% collected from 1978
demnew_NO = 'gimpdem_150m_dem.mat';
demold_NO = 'olddem_1978.mat';
smb_NO = 'smb_anom_Mean_1979_150.mat';
%terminusR_NO = [3195, 2891]; terminusC_NO = [2484, 3119];
terminusR_NO = [3195, 2896, 2372, 2344, 0, 2181, 3594, 0, 5231]; 
terminusC_NO = [2453, 3069, 3839, 4299, 0, 6710, 8164, 0, 8527];
EndR_NO = [0, 3144, 0, 2751, 0, 2663, 3778, 0, 5657];
EndC_NO = [0, 3233, 0, 4318, 0, 6507, 7580, 0, 8507];

NE = [45311, 45310, 70318];
% 45310, 70318 should be in CE but the old dem is collected from 1987 and
% 1985
demnew_NE = 'gimpdem_150m_dem.mat';
demold_NE = 'olddem_1987.mat'; %and 1985
smb_NE = 'smb_anom_Mean_1987_150.mat';
terminusR_NE = [7268, 0]; 
terminusC_NE = [9048, 0];
EndR_NE = [7247, 0];
EndC_NE = [8665, 0];

CE = [44893, 45484, 13995, 20527, 38464];
demnew_CE = 'gimpdem_150m_dem.mat';
demold_CE = 'olddem_1981.mat'; 
smb_CE = 'smb_anom_Mean_1981_150.mat';

SE = [26665, 4677, 23816, 39056, 63077, 41360, 29970, 11855, 25493, 11855];
demnew_CE = 'gimpdem_150m_dem.mat';
%demold_CE = 'olddem_1981.mat'; 
smb_CE = 'smb_anom_Mean_1981_150.mat';

NW = [18971, 59220, 18972, 41515, 71925, 361, 50283, 25089, 17146, 6700, 60835, 33151, 31427, 57168];
demnew_NW = 'gimpdem_150m_dem.mat';
demold_NW = 'olddem_1985.mat'; 
smb_NW = 'smb_anom_Mean_1985_150.mat';

CW = [26614, 55774, 45493, 9788, 56467];
demnew_CW  = 'gimpdem_150m_dem.mat';
%demold_CE = 'olddem_1985.mat'; 
smb_CW = 'smb_anom_Mean_1985_150.mat';

SW = [33910, 47673, 20762, 27660, 32993];
demnew_SW = 'gimpdem_150m_dem.mat';
%demold_CE = 'olddem_1985.mat'; 
smb_SW = 'smb_anom_Mean_1985_150.mat';
%% for finding the centerline
cutoff_mini = -500;
cutoff_max = 1000;
useSavedData = 1;
enlarge_basin = 1;

centerline_method = 'find_centerline_twopoints';
step = 1; %5 % for skiping points
maxdelta = 30;%100; % manully control the difference between the two points determining the centerline points
edge_cutoff = 50;%50; % get rid of the small patch points at each end of the row 
interpDistThreshold = nan;%2000.0;
sample_distance = 50; % interpolating, m
distanceThresh = 50000; 
%% for Pe calculation 
% flow law parameters
n = 3;
m = 1;
% for moving average window
moving_average_thickness_factor = 10;
% for moving average
moving_average_type = 'poly_fit';
% for c_Dprime_D
slope_threshold = 0.0055;%0.01;
flotation_threshold = 75.0;
sea_level = 0.0;
flowLaw = 'effective-pressure';