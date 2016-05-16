%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           volume plotting                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

fprintf('data loading...\n');

base_name_Elmer = {'Elmer'};
base_name_BISICLES = {'BISICLES'};
time_name = {'12to22'};
smb_name = {'20thsmb'};
%beta_name = {'celmer_linear','celmer_stepc95','celmer_stepc2011','cbisicles_linear'};
beta_name = {'celmer_extropln','celmer_extropnofrc'};
beta_legend_name = {'linear\_extrap','no\_friction'};
%temp_name = {'95','2011'};
temp_name = {'95'};

% smb_name = {'elevcorr'};
% beta_name = {'celmer_linear'};
% beta_legend_name = {'E\_linear'};
% temp_name = {'95'};
area_name = {'asf','b3'};
% %model_name = {'B','E'};
% volume_name = {'asf','b3'};
% time_name = {'1990_2011'};
% smb_name = {'elecorr'};
% %smb_name = {'meanmonthly'};
% %Bsmb_name = {'control','90to11'};
% %Esmb_name = {'control','nothermal'};
% beta_name = {'linearbeta'};
% %beta_name = {'linearbeta','95beta','2011beta'};

color_name = ['k','b','k','b'];
dataed = 17;
databd = 17;
grey = [0.7,0.7,0.7];
monthinyr = 12.0;
%legend_posi =  [0.2,0.23,0.1,0.1];
legend_posi = 'SouthWest';
%y_limit = {[-24 1],[-16 1]};
B_path_name = '/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/B_results';
E_path_name = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/E_results';
%B_path_name = '/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/';
%E_path_name = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/data4elmer/';
 figure;
counter = 1;
for q = 1:length(area_name)    
   
   
for k = 1:length(smb_name)

   for i = 1:length(temp_name)
       for j = 1:length(beta_name)
            % Elmer plotting
                
               beta_legend_name_input = cell2mat(beta_legend_name(j));
               
               Edata_fwd = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name(i),'_', beta_name(j),'_',smb_name(k),'_',...
               cell2mat(area_name(q)),'mb.txt'])); 
               Edata_fwd(:,1) = Edata_fwd(:,1) - Edata_fwd(1,1);% the change of the volume
                if strcmp(area_name(q),'asf')
                    H1(counter) = plot(1:length(Edata_fwd(:,1)),Edata_fwd(:,1)/1e9,...
                    ['-',color_name(counter)],'LineWidth',2);
                    legend_E(counter) = {['Austfonna: \beta_{',beta_legend_name_input,'}']};
                    hold on;
                else
                    H1(counter) = plot(1:length(Edata_fwd(:,1)),Edata_fwd(:,1)/1e9,...
                    ['-.',color_name(counter)],'LineWidth',2);
                    legend_E(counter) = {['Basin 3: \beta_{',beta_legend_name_input,'}']};
                    hold on;
                end
               %legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}\_T_{',temp_name_input,'}\_SMB_{',smb_name_input,'}']};
               
               counter = counter + 1;

        end
            
                        
    end
          
end
 
     
         ylabel('Volume (km^3)','FontSize',20);
         %xlabel('Time','FontSize',20);
         xlim([0 120]);
         %ylim([-24 1]);
         grid on;
       
         set(gca,'FontSize',20,'xtick',[0,48,96],...
             'xticklabel',['2012Jan';'2016Jan';'2020Jan']);
         grid on;
         
      

%       legend_handle2 = legend(H,'BISICLES: \beta_l_i_n_e_a_r\_SMB_2_0_c','Elmer/Ice: \beta_l_i_n_e_a_r\_SMB_2_0_c',...
%           'BISICLES: \beta_l_i_n_e_a_r\_SMB_a_b_s','Elmer/Ice: \beta_l_i_n_e_a_r\_SMB_a_b_s',...
%           'BISICLES: \beta_l_i_n_e_a_r\_SMB_a_n_o_m','Elmer/Ice: \beta_l_i_n_e_a_r\_SMB_a_n_o_m',...
%           'BISICLES: \beta_9_5\_SMB_a_b_s','Elmer/Ice: \beta_9_5\_SMB_a_b_s',...
%           'BISICLES: \beta_9_5\_SMB_a_n_o_m','Elmer/Ice: \beta_9_5\_SMB_a_n_o_m',...
%           'BISICLES: \beta_2_0_1_1\_SMB_a_b_s','Elmer/Ice: \beta_2_0_1_1\_SMB_a_b_s',...
%           'BISICLES: \beta_2_0_1_1\_SMB_a_n_o_m','Elmer/Ice: \beta_2_0_1_1\_SMB_a_n_o_m',...
       
%           'Elmer/Ice: \beta_p_o_w_e_r\_SMB_a_b_s');

       A=legend(H1,legend_E,...
          'Location',legend_posi);
        
       set(A,'FontSize',13);
     
       flag = flag+1;
         
 end
