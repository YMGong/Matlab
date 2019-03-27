%************************************************************************%
% Calculating thinning rate
% Check and Modify the configurations in constants_greenland first
% 08/02/2019
% yongmei.gong@vub.be
%************************************************************************%
%clear all;
close all;

%% parameterization
fprintf('parameterization...\n')
%load the constants for greenland
config_greenland  
RegionName = 'NE';
Region = eval(RegionName);
%% Data loading 
fprintf('Data loading ...\n');
load(eval(['smb_',RegionName]));
load(eval(['demold_',RegionName]));
load(eval(['demnew_',RegionName]));
gimpdem_gr=basedata;
%% prepare_all_data
%{
    load Basins_150m.mat
     % 1 km Bedmachine
    xpsn1km = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lon'); xpsn1km = rot90(xpsn1km); % the x coordinates
    ypsn1km = ncread('lonlat_epsg3413_1kmBig_dummy_NObnds.nc','lat'); ypsn1km = rot90(ypsn1km); % the y coordinates
    bed1km = ncread('bed.01.cdf','bed'); bed1km = rot90(bed1km); %
    [lat1km,lon1km] = polarstereo_fwd(ypsn1km,xpsn1km,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
    % 150 m Bedmachine
    xpsn = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','lon'); xpsn = rot90(xpsn); % the x coordinates
    ypsn = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','lat'); ypsn = rot90(ypsn); % the x coordinates
    BMbeddem=ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','bed'); BMbeddem = rot90(BMbeddem);
    mask = ncread('BedMachineGreenland_epsg3413_150m_Withlatlon_inv.nc','mask'); mask = rot90(mask); 
            % 0 = ocean, 1 = ice-free land, 2 = grounded ice, 3 = floating ice, 4 = non-Greenland land
    maskB_fjords = mask; maskB_fjords(mask~=1) = NaN; 
    [latB,lonB] = polarstereo_fwd(ypsn,xpsn,EARTHRADIUS,ECCENTRICITY,LAT_TRUE,LON_POSY); 
    deltLat=latB(1,2)-latB(1,1);deltLon=lonB(1,1)-lonB(2,1);
    [latB_extend,lonB_extend] = meshgrid(latB(1,1)-(500*deltLat):deltLat:latB(end,end)+(500*deltLat),...
                                         lonB(1,1)+(500*deltLon):-deltLon:lonB(end,end)-(500*deltLon)); 
    maskB_fjords_extend = NaN(size(latB_extend));
    maskB_fjords_extend (501:end-500,501:end-500)= maskB_fjords;
    mask_extend = NaN(size(latB_extend));
    mask_extend (501:end-500,501:end-500)= mask;
    BMbeddem_extend = NaN(size(latB_extend));
    BMbeddem_extend(501:end-500,501:end-500)= BMbeddem;
    %BMbeddem_extend(isnan(maskB_fjords_extend))=NaN; % include only the fjords
    Basins_extend = NaN(size(latB_extend));
    Basins_extend(501:end-500,501:end-500)= Basins;
    clear Basins BMbeddem maskB_fjords basedata;
    % dem1980s 150 m 
    demold_gr(demold_gr<0 | demold_gr>4000) = nan;
    demold_gr(mask_extend==0|mask_extend==3|mask_extend==4) = nan;
%}
%% calculte thinning rate
if useSavedData
    fprintf('Loading thinning rate...\n');
    fprintf('**** writting to out.log ****\n');
    load(['Rthinning_',...
        cell2mat(regexp(eval(['demold_',RegionName]),'\d*','Match')),...
        '.mat']);
    fprintf('*****************************\n');
else
    fprintf('Calculating thinning rate...\n');
    fprintf('**** writting to out.log ****\n');
    mask_ice = mask_extend;
    mask_ice(mask_extend~=2) = nan;
    [r_total, r_smb, r_dyn] = get_R_dyn(demold_gr, gimpdem_gr, maskB_fjords_extend, mask_ice, Deltasmb);
    %r_total = get_R_dyn(demold_gr, gimpdem_gr, maskB_fjords_extend);
    fprintf('saving thinning rate...\n');
    save(['Rthinning_',...
        cell2mat(regexp(eval(['demold_',RegionName]),'\d*','Match')),...
        '.mat'],...
        'r_total', 'r_smb', 'r_dyn');
    fprintf('*****************************\n');
end 

%% calculate Pe for different basin
    % now identify the real erosion
    % get only the values in the basins 
    % !!!! maybe select basin !!!!
for i = 1%:length(Region) % do it for different basin
    saveDir = ['./',RegionName, '_', num2str(Region(i))];
    fileID = fopen('out.log','w');
    fprintf(fileID, ['Pe\n',...
        'Region: ',RegionName, ' No. ', num2str(Region(i)), '\n'] );
    fprintf(fileID, '*****************************\n');
    diff_basins = nan(size(r_total));
    diff_basins(Basins_extend == Region(i)) = r_total(Basins_extend == Region(i));
    diff_basins(diff_basins < cutoff_mini | diff_basins > cutoff_max)=nan;
    % attach the the part in Basin that does not have data
    if enlarge_basin
        diff_basins(Basins_extend == Region(i) & isnan(diff_basins)) =...
           Basins_extend(Basins_extend == Region(i) & isnan(diff_basins));
        diff_basins(diff_basins == Region(i))= -999;
    end
    % get rid of artificials
    figure, h=imagesc(diff_basins);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(diff_basins));
    title('150m')
    %}
    %% interpolate to 1 km
    %{
    fprintf('Interpolating...\n');
    xs = -720000;  xe = 1080000;
    ys = -3450000; ye = -530000;
    stepD =1000;
    stepV = 1;
    tmp_array = reshape(latB_extend, size(latB_extend,1)*size(latB_extend,2),1);
    tmp_array(:,2) = reshape(lonB_extend, size(lonB_extend,1)*size(lonB_extend,2),1);
    tmp_array(:,3) = reshape(gimpdem_gr, size(gimpdem_gr,1)*size(gimpdem_gr,2),1);
    [surfaceDEM, ~, ~]  = interp2array(xs,xe,ys,ye,tmp_array,stepV,stepD,'natural');
    %surfaceDEM(surfaceDEM>1.0e3) = 1.0e3;
    figure, h=imagesc(surfaceDEM);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(surfaceDEM));
    title('1km')

    tmp_array = reshape(latB_extend, size(latB_extend,1)*size(latB_extend,2),1);
    tmp_array(:,2) = reshape(lonB_extend, size(lonB_extend,1)*size(lonB_extend,2),1);
    tmp_array(:,3) = reshape(diff_basins, size(diff_basins,1)*size(diff_basins,2),1);
    [diff_basins_1km, ~, ~]  = interp2array(xs,xe,ys,ye,tmp_array,stepV,stepD,'natural');
    diff_basins_1km(diff_basins_1km>1.0e3) = 1.0e3;
    figure, h=imagesc(diff_basins_1km);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(diff_basins_1km));
    title('1km')
    save(['diff_basins_NO_',num2str(Region(i)),'_1km.mat'],'diff_basins_1km');
    

    %load diff_basins_NO_56675_1km.mat
    %load smb_anom_Mean_1979.mat
    %load gimpdem_1km_dem.mat
    % get rid of the unconnected pixels 
    
    % 1km
    r_mask = r_total;
    r_mask(isnan(r_total) | r_mask == 1000)=0.0;
    r_mask = bwareaopen(r_mask,300);
    r_mask = imfill(r_mask,'holes');
    figure, imagesc(r_mask);
    %}

    %% calculate Pe
    fprintf(fileID, 'Finding the centerline...\n');

    rowi = eval(['terminusR_',RegionName]);
    rowi = rowi(i);
    coli = eval(['terminusC_',RegionName]);
    coli = coli(i);
    terminus = [int64(lonB_extend(rowi, coli)), int64(latB_extend(rowi, coli))];
    
    if strcmp(centerline_method, 'find_centerline_twopoints')
        rowi = eval(['EndR_',RegionName]);
        rowi = rowi(i);
        coli = eval(['EndC_',RegionName]);
        coli = coli(i);
        PointE = [int64(lonB_extend(rowi, coli)), int64(latB_extend(rowi, coli))];
    else
        PointE = [NaN, NaN];
    end
    clear rowi coli
    %  A point at the ice front
    % terminus = [-1037000, -362000]; % lon, lat of the terminal point;

    fprintf(fileID, 'The parameters: \n');
    fprintf(fileID, [['terminus = ', num2str(terminus), '\n'],...
        ['end point = ', num2str(PointE), '\n'],...
        ['step = ', num2str(step), '\n'],...
        ['maxdelta = ', num2str(maxdelta), '\n'],...
        ['edge cutoff = ', num2str(edge_cutoff), '\n'],...
        ['interpDistThreshold = ', num2str(interpDistThreshold), '\n'],...
        ['sample distance = ', num2str(sample_distance), '\n'],...
        ['centerline method = ', centerline_method],...
        ['Centerline distance threshold = ', distanceThresh]] );
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    % 150m
    r_mask = diff_basins;
    %r_mask(isnan(diff_basins) | r_mask == -999)=0.0; % don't include the whole basin
    r_mask(isnan(diff_basins))=0.0; % include the whole basin
    r_mask = bwareaopen(r_mask,30000);
    r_mask = imfill(r_mask,'holes');
    
    clear centerPoint
    %bedDEM = bed1km;
    bedDEM = BMbeddem_extend;
    bedDEM(r_mask == 0) = nan;
    surfaceDEM = gimpdem_gr;
    surfaceDEM(r_mask == 0) = nan;
    % somehow lon is y, lat is x...
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  %PointE = [3158, 3242];
%  find_centerline_twopoints(r_mask, lonB_extend, latB_extend,...
%      terminus, PointE, maxdelta);

 
    [centerlinex, centerliney,centerline_distance, surface, bed, dynthinning] =...
        find_centerline(lonB_extend, latB_extend,...
        surfaceDEM, bedDEM,r_dyn,...
        interpDistThreshold, sample_distance,r_mask, step, terminus, PointE,...
        maxdelta, edge_cutoff, distanceThresh,...
        centerline_method);
    
    fprintf(fileID, 'Calculating Pe...\n');
    [Peclet, window]= Pe(centerline_distance, 0, bed, surface,...
        'slope_threshold', slope_threshold,...
        'fileID', fileID);
    fprintf(fileID, 'Calculate Pe moving maxima...\n');
    MM_Peclet = movmax(Peclet, [length(Peclet) 0]);
    [~,MM_PecletIndi,~] = unique(MM_Peclet);
    MM_location = nan(size(Peclet));
    MM_location(MM_PecletIndi) = Peclet(MM_PecletIndi);
    test_plot('Peclet number',centerline_distance,...
        'Peclet', Peclet, 'running maximum', MM_location);  
    %% dynamic thinning rate moving average along the centerline
    fprintf(fileID, 'Calculating dynamic thinning rate moving average along the centerline...\n');
    r_dyn_mvavg = moving_average(dynthinning, window);
    r_dyn_accum = (cumsum(r_dyn_mvavg, 'omitnan')/nansum(r_dyn_mvavg))*100;
    test_plot('dynamic thinning', centerline_distance,...
        'moving average', r_dyn_mvavg, 'cumulative', r_dyn_accum);
    %% find the nearest ocean grid on GISM in the fjord
    fprintf(fileID, 'Finding the fjord mouth...\n');
    resl = lonB_extend(1,1) - lonB_extend(2,1);
    % since the centerline has been extrapolated into the fjord
    % we use the last point to be sure
    mouth_coordi(1) = terminus(1) + resl * ceil((int64(centerlinex(30)) - terminus(1))/resl);
    mouth_coordi(2) = terminus(2) + resl * ceil((int64(centerliney(30)) - terminus(2))/resl);
    % check
    [mouth_coordi(1), mouth_coordi(2)] = find((mouth_coordi(1) ==  int64(lonB_extend)) &...
                      (mouth_coordi(2) ==  int64(latB_extend)));
    figure,imagesc(mask_extend);
    hold on;
    plot(mouth_coordi(2), mouth_coordi(1), 'r.', 'MarkerSize',10);
    fprintf(fileID, '*****************************\n'); 
    %}
    %% save outputs
    fprintf(fileID, ['save to ',saveDir]);
    fclose(fileID);
    system(['rm -rf ', saveDir]);
    system(['mkdir ', saveDir]);
    system(['mv out.log ', saveDir]);
    save([saveDir,'/',RegionName, '_', num2str(Region(i)),'.mat'],...
        'dynthinning','r_dyn_mvavg','r_dyn_accum',...
        'MM_Peclet','Peclet',...
        'mouth_coordi','terminus',...
        'window','centerline_distance');
    save_all_figures(saveDir);
end
fprintf('Finished. :) \n');
    %% plotting
    %{
    fprintf('Plotting...\n')
    figure, h=imagesc(r_total);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(r_total));
    title('r_{total}')

    figure, h=imagesc(r_smb);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(r_smb));
    title('r_{smb}')

    figure, h=imagesc(r_dyn);colormap(jet);colorbar;
    set(h, 'AlphaData',~isnan(r_dyn));
    title('r_{dyn}')
    
    figure, h=imagesc(mask_grdS_extend);
                 colormap(jet);hold on;
                 set(h, 'AlphaData',~isnan(mask_grdS_extend));
                 title('mask_grdS_extend')

    figure, h=imagesc(BMdem_extend);
                 colormap(jet);hold on;
                 set(h, 'AlphaData',~isnan(BMdem_extend));
                 title('BMdem_extend')

    figure, h=imagesc(demold_gr);
                 colormap(jet);hold on;
                 set(h, 'AlphaData',~isnan(demold_gr));      

    figure, h=imagesc(gimpdem_gr);
                 colormap(jet);hold on;
                 set(h, 'AlphaData',~isnan(gimpdem_gr));                   

    figure, h=imagesc(diff_dem_gr);
                 colormap(jet);hold on;
                 set(h, 'AlphaData',~isnan(diff_dem_gr));
                 title('diff_demP');

    maskB_fjords_extend(isnan(maskB_fjords_extend))=0.0;             
    figure, h=imagesc(diff_basins);
                 colormap(jet);colorbar;hold on;
                 set(h, 'AlphaData',~isnan(diff_basins));
    %              contour(maskB_fjords_extend, 1);
                 title('diff_basins');

    % figure, h=imagesc(diff_dem_gr);
    %              colormap(jet);hold on;
    %              set(h, 'AlphaData',~isnan(diff_dem_gr));
    %              title('diff_demold_gr')
    %}
    %%
% functions
function data_mvavg = moving_average(data, window)
    data_mvavg = nan(size(data));
    for idx = 1 : length(data)
        if isnan(window(idx))
            continue
        end
             % Find elements within window
          idx_start = int64(max(0     , floor(idx - window(idx)/2)) );
          idx_stop  = int64(min(length(data), floor(idx + window(idx)/2)) );
          data_window   = data(idx_start:idx_stop);
          if any(~isnan(data_window))
            data_mvavg(idx) =  nanmean(data_window);
          else
            data_mvavg(idx) =  nan;
          end
    end
    data_mvavg(isnan(data)) = nan;
end
function save_all_figures(saveDir)
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
      FigHandle = FigList(iFig);
      FigName   = get(FigHandle, 'Number');
      savefig(FigHandle, fullfile(saveDir, ['/figure_',num2str(FigName), '.fig']));    
    end
end
