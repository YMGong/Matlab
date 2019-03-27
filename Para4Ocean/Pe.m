 function [Peclet, window] = Pe(centerline_distance, terminus, bed, surface, varargin)
%% Get the parameters
format long
config_greenland

    for i = 1:2:length(varargin)
            switch varargin{i}
                case 'm'
                      m = varargin{i+1};
                case 'n'
                      n = varargin{i+1};
                case 'moving_average_thickness_factor'
                    moving_average_thickness_factor = varargin{i+1};
                case 'moving_average_type'
                    moving_average_type = varargin{i+1};
                case 'slope_threshold'
                    slope_threshold = varargin{i+1};
                case 'flotation_threshold'
                    flotation_threshold = varargin{i+1};
                case 'sea_level'
                    sea_level = varargin{i+1};
                case 'flowLaw'
                    flowLaw = varargin{i+1};
                case 'fileID'
                    fileID = varargin{i+1};
                otherwise
                      error(['Unexpected option: ' varargin{i}])
            end
    end

fprintf(fileID, 'The parameters: \n');
fprintf(fileID, [['m = ', num2str(m), '\n'],...
        ['n = ', num2str(n), '\n'],...
        ['moving_average_thickness_factor = ', num2str(moving_average_thickness_factor), '\n'],...
        ['slope_threshold = ', num2str(slope_threshold), '\n'],...
        ['flotation_threshold = ', num2str(flotation_threshold), '\n'],...
        ['sea_level = ', num2str(sea_level), '\n'],...
        ['flowLaw = ', flowLaw, '\n']]);
%% Data preparation 
fprintf(fileID, 'Smoothing surface and bed...\n');
x = centerline_distance;
% Length of perturbation (l)
% l is just the distance between two centerline points, so that the slope
% alpha = H/L can be estimated
l = x;
l(centerline_distance < terminus) = nan;  % nan out everything downglacier of terminus
l = l - nanmin(l);                        % shift l such that zero is at terminus

% ice thickness
h  = surface - bed;
% water height;
hw = sea_level - bed;

window = moving_average_window_by_thickness(centerline_distance,...
                                            bed,...
                                            surface,...
                                            moving_average_thickness_factor);                                        
surface = moving_fit(centerline_distance, surface, window);
bed = moving_fit(centerline_distance, bed, window);
test_plot('Profile along centerline after fitting',centerline_distance,...
    'surface', surface,...
    'bed', bed);

%% Calculate moving average window
moving_average_window = moving_average_window_by_thickness(centerline_distance,...
                                                              bed,...
                                                              surface,...
                                                              moving_average_thickness_factor);
%%  Moving fit
[hfit, dhfit, slopefit, dslopefit, hwfit, dhwfit] = Pe_moving_average(x, surface, h,...
      hw, moving_average_window);
  
test_plot('Ice thickness along centerline',centerline_distance,...
    'fitted h', hfit,...
    'original h', h);  
% Set water depth to 0.0 where bed goes above sea level
   nanidx = find(isnan(hwfit));
   
   hwfit(nanidx) = -9999.0;
   dhwfit(nanidx) = -9999.0;
   hwfit(hwfit < 0) = 0.0;
   dhwfit(hwfit < 0) = 0.0;
   hwfit(nanidx) = nan;
   dhwfit(nanidx) = nan;

   
%% Kinematic wave coefficients
% 'hard-bed'
%[c, D, Dprime]=c_Dprime_D(hfit, dhfit, slopefit, dslopefit, flowLaw, slope_threshold, n);

% 'effective-pressure'
[c, D, Dprime]=c_Dprime_D(hfit, dhfit, slopefit, dslopefit, flowLaw, slope_threshold,...
                            'hw', hwfit, 'dhw', dhwfit, 'flotation_threshold', flotation_threshold,...
                            'centerline_distance', centerline_distance);

%% Peclet number
Peclet = ( (c - Dprime) ./ D ) .* l;
%% Remove values before/after discontinuities until slope break
% maybe it is not nessary ...

end
function  window = moving_average_window_by_thickness(centerline_distance, bed, surface, thickness_factor)
   
   % thickness_factor is 10 in the paper 
   % bed and surface is the 'true' value along the centerline %
   h = surface - bed;
   window = nan(size(h));
   mvavgwindow = floor(thickness_factor * h);
   for idx = 1: length(h)
      window(idx) = mvavgwindow(idx) / nanmean(diff(centerline_distance));   
      % only average for the value at the centerline (have not taken into account the value outside the centerline) 
      % Round n to odd number
      window(idx) = ceil(window(idx) / 2.) * 2 - 1;
   end
   valididx = find(~isnan(window));
   xvalid = centerline_distance(~isnan(window));
   %windowvalid = window(~isnan(window));
end
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
function [data_mvavg, delta_mvavg] = moving_fit(x, data, window)
data_mvavg = nan(size(data));

    for idx = 1 : length(data)
        if isnan(window(idx))
            continue
        end
             % Find elements within window
        idx_start = int64(max(0     , floor(idx - window(idx)/2)) );
        idx_stop  = int64(min(length(data), floor(idx + window(idx)/2)) );
        
        if idx_start~=0 && idx_stop~=0
            x_window   = x(idx_start:idx_stop);
            data_window   = data(idx_start:idx_stop);
            % Find non-nan entries
            valididx = find(~isnan(data_window));
            x_valid = x_window(valididx);
            data_valid = data_window(valididx);
            % fit into a polynomial function
            n = length(x_valid);
            if n > 2
               %[p1,S1] = polyfit(x_valid, data_valid, 1);
                H = [(x_valid.'), ones(size(x_valid)).'];
                V = inv( transpose(H) * H );
                p = V * transpose(H) * (data_valid.');
                
                data_mvavg(idx) = polyval(p, x(idx));
                
            end
            
        end
    end
    
    data_mvavg(isnan(data)) = nan;
end
function [data_fit, data_fiterr, dy_fit, dy_fiterr] =...
    fit_first_order(x, data, varargin)

% shift coordinates
   x_shift = mean(x);
   x = x - x_shift;
   data_shift = mean(data);
   data = data - data_shift;
   
   n = length(data);
       if n > 3
           
           %[p,S] = polyfit(x, data, data_err);
           H = [(x.'), ones(size(x)).'];
           
           if exist('data_err', 'var')
               W = diag(1. ./ data_err.^2 );
               V = inv( transpose(H) * W * H );
               p = V * transpose(H) * W * (data.');
           else
               V = inv( transpose(H) * H );
               p = V * transpose(H) * (data.');
           end   
 
           if strcmp(varargin{1}, 'x_fit')
               x_fit = varargin{2} - x_shift;
           else
               x_fit = x;
           end
           
           data_fit = polyval(p, x_fit);
           data_fit = data_fit + data_shift;
           data_fiterr   = sqrt(V(2,2)); % TBD: This should probably be the RSS of the post-fit residuals
           %hw_fiterr  = np.sqrt(V(2,2)); % TBD: This should probably be the RSS of the post-fit residuals
           dy_fit     = p(1);
           dy_fiterr  = sqrt(V(1,1)); % TBD: What should this be?

       end
            
end
function [data_fit, data_fiterr, ddata_fit, ddata_fiterr, dddata_fit, dddata_fiterr] =...
    fit_second_order(x, data, varargin)
% shift coordinates
   x_shift = mean(x);
   x = x - x_shift;
   data_shift = mean(data);
   data = data - data_shift;
   
   n = length(data);
   if n > 4
       H = [(x.').^2, (x.'), ones(size(x)).'];
       
       if exist('data_err', 'var')
           W = diag(1. ./ data_err.^2 );
           V = inv( transpose(H) * W * H );
           p = V * transpose(H) * W * (data.');
       else
           V = inv( transpose(H) * H );
           p = V * transpose(H) * (data.');
           
       end    
       
       if strcmp(varargin{1}, 'x_fit')
         x_fit = varargin{2} - x_shift;
       else
         x_fit = x;
       end
       data_fit = polyval(p, x_fit);
       data_fit = data_fit + data_shift;   
       data_fiterr      = sqrt(V(3,3));      % TBD: Placeholder - this needs to be propagated through second order polynomial
                                          % => Var[sfit] = [1 x x^2] * V * [1 x x^2]'
       ddata_fit        = polyval(p(1:2), x_fit);
%        tmp = polyfit(x,data,1);
%        ddata_fit        = tmp(:, 1);
       ddata_fiterr     = sqrt(V(2,2)); % TBD: Placeholder - this needs to be propagated through second order polynomial
       dddata_fit       = p(1);
       dddata_fiterr    = sqrt(V(1,1));% TBD: Placeholder - this needs to be propagated through second order polynomial
   end
   
end
function [hfit, dhfit, slopefit, dslopefit, hwfit, dhwfit] = ...
    Pe_moving_average(x, surface, thickness, water_depth, moving_average_window)
    %% Initialize
   sfit         = nan(size(x));
   sfiterr      = nan(size(x));
   hfit         = nan(size(x));
   hfiterr      = nan(size(x));
   hwfit        = nan(size(x));
   hwfiterr     = nan(size(x));
   dhfit        = nan(size(x));
   dhfiterr     = nan(size(x));
   dhwfit       = nan(size(x));
   dhwfiterr    = nan(size(x));
   slopefit     = nan(size(x));
   slopefiterr  = nan(size(x));
   dslopefit    = nan(size(x));
   dslopefiterr = nan(size(x));
   
   s     = surface;
%    serr  = unp.std_devs(surface)
   h     = thickness;
%    herr  = unp.std_devs(thickness)
   hw    = water_depth;
%    hwerr = unp.std_devs(water_depth)
%% use fit function by default
   % OPTION fit polynomials but NOT at the end points, where window extends beyond valid data
   % This is also called a savitzky_golay filter (there's a python function to do this)  
   firstValid = find(~isnan(surface), 1, 'first');
   lastValid = find(~isnan(surface), 1, 'last');
    for idx = firstValid : lastValid
        % Find elements within window
        idx_start = int64(max(0     , floor(idx - moving_average_window(idx)/2)) );
        idx_stop  = int64(min(length(h), floor(idx + moving_average_window(idx)/2)) );
        if idx_start~=0 && idx_stop~=0
            valididx  = find(~isnan(s(idx_start:idx_stop)));
            % Find valid entries within the window
            xwindowvalid  = x(idx_start:idx_stop);
            xwindowvalid = xwindowvalid(valididx);
            swindowvalid  = s(idx_start:idx_stop);
            swindowvalid = swindowvalid(valididx);
            hwindowvalid  = h(idx_start:idx_stop);
            hwindowvalid = hwindowvalid(valididx);
            hwwindowvalid  = hw(idx_start:idx_stop);
            hwwindowvalid = hwwindowvalid(valididx);
          
            if idx_start >= firstValid && idx_stop <= lastValid
                
               [sfit(idx), sfiterr(idx), slopefit(idx), slopefiterr(idx), dslopefit(idx), dslopefiterr(idx)] =...
                    fit_second_order(xwindowvalid, swindowvalid, 'x_fit',  x(idx));
                
               [hfit(idx),  hfiterr(idx),  dhfit(idx),     dhfiterr(idx)]  =...
                   fit_first_order( xwindowvalid, hwindowvalid, 'x_fit',  x(idx));
               
               [hwfit(idx), hwfiterr(idx), dhwfit(idx),    dhwfiterr(idx)] =...
                   fit_first_order( xwindowvalid, hwwindowvalid, 'x_fit',  x(idx));
               
               
            end
        end
    end
    invalididx = find(isnan(surface));
    sfit(invalididx) = nan;      sfiterr(invalididx) = nan;
    slopefit(invalididx) = nan;  slopefiterr(invalididx) = nan;
    dslopefit(invalididx) = nan; dslopefiterr(invalididx) = nan;
    hfit(invalididx) = nan;      hfiterr(invalididx) = nan;
    dhfit(invalididx) = nan;     dhfiterr(invalididx) = nan;
    hwfit(invalididx) = nan;      hwfiterr(invalididx) = nan;
    dhwfit(invalididx) = nan;     dhwfiterr(invalididx) = nan;
 
    
end

function [c, D, Dprime]=c_Dprime_D(h,dh,slope,dslope, flowLaw, slope_threshold, varargin)
   %Note: "Kb" has been factored out (see Cuffey and Paterson, 2010)
   config_greenland 
   c = nan(size(slope));
   D = nan(size(slope));
   Dprime = nan(size(slope));
   % default
   n = 3.0;
   m = 1.0;
   
   %% calculate according to sliding law
   if strcmp(flowLaw, 'hard-bed') 
          
      for i = 1:2:length(varargin)
        switch lower(varargin{i})
              case 'm'
                  m = varargin{i+1};
%               case 'n'
%                   n = varargin{i+1};
              otherwise
                  error(['Unexpected option: ' varargin{i}])
        end
      end
      
	  % Invalidate where slopes are too small
      slope(slope <= slope_threshold) = nan;
      validIdx = find(~isnan(slope));
      % Note: "Kb" has been factored out (see Cuffey and Paterson, 2010)
      c(validIdx) = (m+1) * (h(validIdx) * slope(validIdx))^m;
      D(validIdx) = m * h(validIdx)^(m+1) .* slope(validIdx)^(m-1);
      Dprime(validIdx) = m * (m+1) * h(validIdx) ^ m .* dh(validIdx) .* slope(validIdx) ^ (m-1) ...
                         + m * (m-1) * h(validIdx)^(m+1) .* slope(validIdx)^(m-2) .* dslope(validIdx);
   elseif strcmp(flowLaw, 'effective-pressure')
      flotation_threshold = 75.0;
      
      for i = 1:2:length(varargin)
        switch lower(varargin{i})
            case 'm'
                m = varargin{i+1};
            case 'n'
                n = varargin{i+1};
            case 'hw'
                hw = varargin{i+1};
            case 'dhw'
                dhw = varargin{i+1};
            case 'flotation_threshold'
                flotation_threshold = varargin{i+1};
            case 'centerline_distance'
                centerline_distance = varargin{i+1};
            otherwise
                error(['Unexpected option: ' varargin{i}])
        end
      end
      
      phatw = (rhow/rhoi) * hw;
      dphatw = (rhow/rhoi) * dhw;
      test_plot('slope along centerline',centerline_distance,...
    'slope', slope);  
          validIdx = find(((h - phatw) > flotation_threshold) & (slope > slope_threshold));
          %invalidIdx = np.logical_not(validIdx);

          % Note: "kappa" has been factored out (see Pfeffer, 2007)
          peff = h - phatw;
          c(validIdx) = ( ((h(validIdx).* slope(validIdx)).^ n) ./ (peff(validIdx) .^ m) ) .*...
                        ( 1 + ((n-m)*h(validIdx) - n*phatw(validIdx)) ./ peff(validIdx) );
                    
          D(validIdx) = n * ( (slope(validIdx) .^ (n-1) .* h(validIdx) .^ (n+1)) ./ peff(validIdx) .^ m );
          
          Dprime(validIdx) = ( n*(n-1)*slope(validIdx) .^ (n-2) .* dslope(validIdx) .* h(validIdx) .^ (n+1) ) / ( peff(validIdx) .^ m ) + ...
                             ( n*(n+1)*slope(validIdx) .^ (n-1) .* dh(validIdx) .* h(validIdx) .^ n ) ./ ( peff(validIdx) .^ m ) - ...
                             ( n*m*(dh(validIdx)-dphatw(validIdx)).* slope(validIdx) .^ (n-1) .* h(validIdx) .^ (n+1) ) ./ ( peff(validIdx) .^ (m+1) );
      
   else
       error('Only hard-bed and effective-pressure flow laws can be chosen.\n')
   end
   
end

% function add_uncertainty
% end
