function [r_total, r_smb, r_dyn] = get_R_dyn(demold, demnew, mask_rock, mask_ice, r_smb, varargin)
%function r_total = get_R_dyn(demold, demnew, mask_rock, varargin)
    fprintf('Calculating gradient...\n')
    % construction intervals for the mask
    if ~isempty(varargin)
        for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'grd_intervl'
                          grd_intervl = varargin{i+1};
                    otherwise
                          error(['Unexpected option: ' varargin{i}])
                end
        end
    else
        grd_intervl = [0 2 4 8 15 30 50]; % the default

    end
    
    [dsx, dsy] = gradient(demold); % according to the old dem 
    mask_grdS = sqrt(dsx.^2+dsy.^2);

    for i = 1:length(grd_intervl)
        if i+1>length(grd_intervl)
            mask_grdS(mask_grdS>=grd_intervl(i))=i;
        else
            mask_grdS(mask_grdS>=grd_intervl(i) &...
                mask_grdS<grd_intervl(i+1))=i;
        end
    end
    mask_grdS_rock = mask_grdS;
    mask_grdS_rock(isnan(mask_rock)) = nan;
    fprintf('Please Check if the points are evenly distributed in each intervals!\n');
    figure,h=imagesc(mask_grdS); colormap(jet),colorbar; set(h, 'AlphaData',~isnan(mask_grdS));
    clear dsx dsy;
    %%
    fprintf('Correcting...\n');
    % To remove the bias, vertical differences between GIMP DEM and 
    % each 1985 DEM were found over bedrock, as defined by the BedMechine surface mask, for slopes less than 20?
    % and the mean vertical offset was removed from the GIMP DEM.
    demnew(isnan(demold)) = nan;
    %BMdem_extend(isnan(demold)|isnan(gimpdem_gr))=nan;

    fprintf('Eliminate Outliers Using Interquartile Range then calculate the mean...\n')
    demdiff = demnew - demold;
    
    for i = 1:max(max(mask_grdS_rock))
        demdiff_tmp = demdiff;
        demdiff_tmp(isnan(mask_grdS_rock)| mask_grdS_rock~=i)=NaN;
        %demdiff = demdiff_tmp(~isnan(demdiff_tmp));
        if ~isnan(max(max(demdiff_tmp)))
            % calculate the adjustment over the rocks
            demdiffP = demdiff_tmp(~isnan(demdiff_tmp) & demdiff_tmp>=0);
            demdiffN = demdiff_tmp(~isnan(demdiff_tmp) & demdiff_tmp<0);
            mean_diffP = mean_elim(demdiffP); % eliminate outliers for the positive end
            mean_diffN = mean_elim(demdiffN); % eliminate outliers for the positive end

            % adjust the vertical difference of both rocks and ice
            % according to the nagtive and positive adjustment
            demdiff((demdiff>=0) & mask_grdS == i) =...
                demdiff(demdiff>=0 & mask_grdS == i) + mean_diffP;

            demdiff(demdiff<0 & mask_grdS == i) =...
                demdiff(demdiff<0 & mask_grdS == i) + mean_diffN; 
        end
    end
    figure,h=imagesc(demdiff); colormap(jet),colorbar; set(h, 'AlphaData',~isnan(demdiff));
    
    
    %% r_dyn
    r_total = demdiff;
    r_total(isnan(mask_ice))=nan;
    r_dyn = r_total - r_smb;
    %% 
end
function avg = mean_elim(array)
    array = sort(array);
    avg = mean(array);
    trd = 1.5*iqr(array);

    while (avg-array(1))>trd
        array=array(2:end);
        avg = mean(array);
        trd = 1.5*iqr(array);
        while (array(end)-avg)>trd
            array=array(1:end-1); 
            avg = mean(array);
            trd = 1.5*iqr(array);
        end   
    end
end
 