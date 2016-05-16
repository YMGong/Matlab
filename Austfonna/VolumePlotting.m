%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  volume plotting                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all;

fprintf('data loading...\n');

base_name_Elmer = {'Elmer'};
base_name_BISICLES = {'BISICLES'};
time_name = {'95to11'};
smb_name = {'elevcorr','noelevcorr'};
%beta_name = {'celmer_linear','celmer_stepc95','celmer_stepc2011','cbisicles_linear'};
beta_name = {'celmer_linear','celmer_stepc95','celmer_stepc2011'};
%beta_legend_name = {'E\_linear','E\_step95','E\_step2011','B\_linear'};
beta_legend_name = {'linear','step1995','step2011'};
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

color_name = ['k','b','m','r','c','y'];
dataed = 17;
databd = 17;
grey = [0.7,0.7,0.7];
monthinyr = 12.0;
%legend_posi =  [0.2,0.23,0.1,0.1];
legend_posi = 'SouthWest';
%y_limit = {[-24 1],[-16 1]};
B_path_name = 'J:\ASF_data\B_results';
E_path_name = 'J:\ASF_data\E_results';
%B_path_name = '/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/B_results';
%E_path_name = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/E_results';
%B_path_name = '/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/';
%E_path_name = '/wrk/ygong/ASFforwardFiles/Invers/postprocessing/data4elmer/';


for q = 1:length(area_name)    
    figure;
    counter = 1;flag = 1;
    counter2 = 1;
    
        Bfile_name = cell2mat([B_path_name,'/',base_name_BISICLES,'_',...
            time_name,'_temp',cell2mat(temp_name(1)),'_',...
            cell2mat(beta_name(1)),'_20thsmb_',cell2mat(area_name(q))]);
        Efile_name = cell2mat([E_path_name,'/',base_name_Elmer,'_',...
            time_name,'_temp',cell2mat(temp_name(1)),'_',...
            cell2mat(beta_name(1)),'_20thsmb_',cell2mat(area_name(q))]);

        
        Bdata =load([Bfile_name,'mb.tab']);
        [tmp, index] = sort(Bdata(:,1),1);
        Bdata(:,1)=tmp;
        Bdata(:,2)=Bdata(index,2);
        Bdata(:,3)=Bdata(index,3);
        Bdata_change = ((Bdata(end,1)-Bdata(1,1))/Bdata(1:1))*100;
        Bdata(:,2) = Bdata(:,2) - Bdata(1,2); % the change of the volume
        Edata =load([Efile_name,'mb.txt']);
        Edata_change = (Edata(end,1)-Edata(1,1))/Edata(1:1)*100;
        Edata(:,1) = Edata(:,1) - Edata(1,1); % the change of the volume
           
        H1(counter)=plot(1:length(Edata(:,1)),Edata(:,1)/1e9,'-','Color',grey,'LineWidth',2);
        %legend_E(counter) = {['Elmer/Ice: \beta_{',cell2mat(beta_legend_name(1)),'}\_T_{',cell2mat(temp_name(1)),'}\_SMB_{20th}']};
        legend_E(counter) = {['Elmer/Ice: \beta_{',cell2mat(beta_legend_name(1)),'}','\_SMB_{90s}']};
        counter = counter + 1;
        hold on;
        
%         %H2(counter2)=plot(0:1:monthinyr*(databd)-1,Bdata(:,2)/1e9,'-.','Color',grey,'LineWidth',2);
%         H2(counter2)=plot(Bdata(:,1)*12.0,Bdata(:,2)/1e9,'-.','Color',grey,'LineWidth',2);
%         %legend_B(counter2) = {['BISICLES: \beta_{',cell2mat(beta_legend_name(1)),'}\_T_{',cell2mat(temp_name(1)),'}\_SMB_{20th}']};
%         legend_B(counter2) = {['BISICLES: \beta_{',cell2mat(beta_legend_name(1)),'}','\_SMB_{90s}']};
%         counter2 = counter2 + 1;
        
        hold on;
    
 for k = 1:length(smb_name)

    for i = 1:length(temp_name)
        for j = 1:length(beta_name)
            % Elmer plotting
            if strcmp(smb_name(k),'elevcorr') && strcmp(temp_name(i),'95')
               smb_name_input = 'elevcorr';
               temp_name_input = '95';
               beta_name_input = cell2mat(beta_name(j));
               beta_legend_name_input = cell2mat(beta_legend_name(j));
               
               Edata_fwd = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name_input,'_', beta_name_input,'_',smb_name_input,'_',...
               cell2mat(area_name(q)),'mb.txt'])); 
               Edata_fwd(:,1) = Edata_fwd(:,1) - Edata_fwd(1,1);% the change of the volume
    
               H1(counter) = plot(1:length(Edata(:,1)),Edata_fwd(:,1)/1e9,...
              ['-',color_name(counter-1)],'LineWidth',2);
               hold on;
               
               %legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}\_T_{',temp_name_input,'}\_SMB_{',smb_name_input,'}']};
               legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}','\_SMB_{',smb_name_input,'}']};
               counter = counter + 1;
            elseif strcmp(smb_name(k),'elevcorr') && strcmp(temp_name(i),'2011') && strcmp(beta_name(j),'celmer_linear')
               smb_name_input = 'elevcorr';
               temp_name_input = '2011';
               beta_name_input = 'celmer_linear'; 
               beta_legend_name_input = cell2mat(beta_legend_name(1));
               
               Edata_fwd = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name_input,'_', beta_name_input,'_',smb_name_input,'_',...
               cell2mat(area_name(q)),'mb.txt'])); 
               Edata_fwd(:,1) = Edata_fwd(:,1) - Edata_fwd(1,1);% the change of the volume
    
               H1(counter) = plot(1:length(Edata(:,1)),Edata_fwd(:,1)/1e9,...
              ['-',color_name(counter-1)],'LineWidth',2);
               hold on;
               
               %legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}\_T_{',temp_name_input,'}\_SMB_{',smb_name_input,'}']};
               legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}','\_SMB_{',smb_name_input,'}']};
               counter = counter + 1;
            elseif strcmp(smb_name(k),'noelevcorr') && strcmp(temp_name(i),'95') && strcmp(beta_name(j),'celmer_linear') 
               beta_name_input = 'celmer_linear'; 
               temp_name_input = '95';
               smb_name_input = 'noelevcorr';
               beta_legend_name_input = cell2mat(beta_legend_name(1));
               
               Edata_fwd = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name_input,'_', beta_name_input,'_',smb_name_input,'_',...
               cell2mat(area_name(q)),'mb.txt'])); 
               Edata_fwd(:,1) = Edata_fwd(:,1) - Edata_fwd(1,1);% the change of the volume
    
               H1(counter) = plot(1:length(Edata(:,1)),Edata_fwd(:,1)/1e9,...
             ['-',color_name(counter-1)], 'LineWidth',2);
               hold on;
               
               %legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}\_T_{',temp_name_input,'}\_SMB_{',smb_name_input,'}']};
               legend_E(counter) = {['Elmer/Ice: \beta_{',beta_legend_name_input,'}','\_SMB_{',smb_name_input,'}']};
               counter = counter + 1;
            end
            
            %BISICLES plotting
            if strcmp(smb_name(k),'elevcorr') && strcmp(beta_name(j),'celmer_linear')
                smb_name_Binput = 'elevcorr';
                beta_name_Binput = 'celmer_linear';
                temp_name_Binput = cell2mat(temp_name(i));
                beta_legend_name_input = cell2mat(beta_legend_name(1));
                
                Bdata_fwd = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
                temp_name_Binput,'_', beta_name_Binput,'_',smb_name_Binput,'_',...
                cell2mat(area_name(q)),'mb.tab'])); 
                Bdata_fwd(:,2) = Bdata_fwd(:,2) - Bdata_fwd(1,2);% the change of the volume
               
                H2(counter2) = plot(Bdata(:,1)*12.0,Bdata_fwd(:,2)/1e9,...
                ['-.',color_name(counter2)], 'LineWidth',2); 
               %['-.',color_name(counter2-1)], 'LineWidth',2);
                hold on;
                
                %legend_B(counter2) = {['BISICLES: \beta_{',beta_legend_name_input,'}\_T_{',temp_name_Binput,'}\_SMB_{',smb_name_Binput,'}']};
                legend_B(counter2) = {['BISICLES: \beta_{',beta_legend_name_input,'}','\_SMB_{',smb_name_Binput,'}']};
                counter2 = counter2 + 1;
            end
              
        end
          
    end
 end
      if strcmp(area_name(q),'asf')
         ylabel('Volume (km^3)','FontSize',15);
         %xlabel('Time','FontSize',20);
         xlim([0 204]);
         ylim([-24 1]);
         grid on;
       
         set(gca,'FontSize',15,'xtick',[0,48,96,144,192],...
             'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan']);
         title('Volume change of ASF');
         %text(20,-5,['BISICLES:mass change < ',num2str(Bdata_change),'%'],'FontSize',15);
         %text(20,20,['Elmer/Ice:mass change < ',num2str(Edata_change),'%'],'FontSize',15);
         grid on;
         %lineH = line([monthinyr*(databd-1) monthinyr*(databd-1)],[-30 40]);
         %text(194,37,'Relaxation','FontSize',15);text(235,37,'Present Day','FontSize',15);
         %set(gca,'XTickLabel',[0 5 10 15 20 1994 1999 2004 2009 2014]);
      
      else
        ylabel('Volume (km^3)','FontSize',15);
        %xlabel('Time','FontSize',20);
        xlim([0 204]);
        ylim([-16 1]);
        set(gca,'FontSize',15,'xtick',[0,48,96,144,192],...
            'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan']);
        title('Volume change of Basin 3');
        %text(6,6,['BISICLES: mass change <',num2str(Bdata_change),'%'],'FontSize',15);
        %text(6,-2,['Elmer/Ice: mass change <',num2str(Edata_change),'%'],'FontSize',15);
        grid on;
        %lineH = line([monthinyr*(databd-1) monthinyr*(databd-1)],[-10 10]);
        %text(194,9,'Relaxation','FontSize',15);text(235,9,'Present Day','FontSize',15);
        %set(gca,'XTickLabel',[0 10 20 1994 2004 2014])
      end
     if q==1
      %H = gca;  
      H = [rot90(H1,-1);rot90(H2,-1)];
%       legend_handle2 = legend(H,'BISICLES: \beta_l_i_n_e_a_r\_SMB_2_0_c','Elmer/Ice: \beta_l_i_n_e_a_r\_SMB_2_0_c',...
%           'BISICLES: \beta_l_i_n_e_a_r\_SMB_a_b_s','Elmer/Ice: \beta_l_i_n_e_a_r\_SMB_a_b_s',...
%           'BISICLES: \beta_l_i_n_e_a_r\_SMB_a_n_o_m','Elmer/Ice: \beta_l_i_n_e_a_r\_SMB_a_n_o_m',...
%           'BISICLES: \beta_9_5\_SMB_a_b_s','Elmer/Ice: \beta_9_5\_SMB_a_b_s',...
%           'BISICLES: \beta_9_5\_SMB_a_n_o_m','Elmer/Ice: \beta_9_5\_SMB_a_n_o_m',...
%           'BISICLES: \beta_2_0_1_1\_SMB_a_b_s','Elmer/Ice: \beta_2_0_1_1\_SMB_a_b_s',...
%           'BISICLES: \beta_2_0_1_1\_SMB_a_n_o_m','Elmer/Ice: \beta_2_0_1_1\_SMB_a_n_o_m',...
       
%           'Elmer/Ice: \beta_p_o_w_e_r\_SMB_a_b_s');

       A=legend(H,cell2mat(legend_E(1)),cell2mat(legend_E(2)),...
          cell2mat(legend_E(3)),cell2mat(legend_E(4)),...
          cell2mat(legend_E(5)),...
          cell2mat(legend_B(1)),...
          'Location',legend_posi);
        
       set(A,'FontSize',10);
     end
       flag = flag+1;
         
       %set(lineH,'LineWidth',2,'LineStyle','-.','Color','k');
 end
      

      
% print(gcf,'-dtiff','VolumeChange_BISICLES.tiff');
%saveas(gcf,['.\pics\','volunm_change_',cell2mat(volume_name(k)),'.tif'],'tif');


%set(lineH,'LineWidth',2,'LineStyle','-.','Color','k');
if 0
hold on;
fprintf('plotting...\n');
data=load('/wrk/ygong/BISICLES/ASF_simulation/PostProcessing/test_b3mb.tab');
plot(0:length(data)-1,(data(:,2)-data(1,2))/1e9,'.r');

figure;
Edata_500yr = load([E_path_name,'/',base_name_Elmer,...
              'linearbeta_ebet_',cell2mat(time_namt(2)),Efile_name,'.txt']);
      Edata_500yr(:,1) = Edata_500yr(:,1) - Edata_500yr(1,1); 
      plot(dataed+1:dataed+length(Edata_500yr),Edata_500yr(:,1)/1e9,...
              '-','Color',grey,'LineWidth',2),hold on;   
end   

fprintf('plotting...\n');

 ylimit={[18.58 18.64],[18.36 18.42]};
 xlimit=[0 204];
 Edata_fwd_asf = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
               cell2mat(area_name(1)),'mb.txt']));
 Edata_fwd_b3 = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
               cell2mat(area_name(2)),'mb.txt']));
           
 Bdata_fwd_asf = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
                temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
                cell2mat(area_name(1)),'mb.tab']));  
 [tmp, index] = sort(Bdata_fwd_asf(:,1),1);
        Bdata_fwd_asf(:,1)=tmp;
        Bdata_fwd_asf(:,2)=Bdata_fwd_asf(index,2);
        Bdata_fwd_asf(:,3)=Bdata_fwd_asf(index,3);           
 Bdata_fwd_b3 = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
                temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
                cell2mat(area_name(2)),'mb.tab']));  
 [tmp, index] = sort(Bdata_fwd_b3(:,1),1);
        Bdata_fwd_b3(:,1)=tmp;
        Bdata_fwd_b3(:,2)=Bdata_fwd_b3(index,2);
        Bdata_fwd_b3(:,3)=Bdata_fwd_b3(index,3);           
            
figure,plot(0:1:monthinyr*(dataed)-1,abs((Edata_fwd_b3(:,1)-Edata_fwd_b3(1,1))./abs(Edata_fwd_asf(:,1)-Edata_fwd_asf(1,1)))*100,...
              '-k','LineWidth',2);hold on;
plot(Bdata_fwd_b3(:,1)*12.0,abs((Bdata_fwd_b3(:,2)-Bdata_fwd_b3(1,2))./abs(Bdata_fwd_asf(:,2)-Bdata_fwd_asf(1,2)))*100,...
              '-r','LineWidth',2);hold on;grid on;
xlim(xlimit);
% set(gca,'xtick',[0,48,96,144,192],...
%              'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],'FontSize',15);
         
legend('Elmer','BISICLES');title('\Delta(V_{b3}/V_{ASF})x100%');

figure,[AX,Hplot1,Hplot2]=plotyy(0:1:monthinyr*(dataed)-1,(Edata_fwd_b3(:,1)./Edata_fwd_asf(:,1))*100,...
    Bdata_fwd_b3(:,1)*12.0,(Bdata_fwd_b3(:,2)./Bdata_fwd_asf(:,2))*100);
tmp=cell2mat(ylimit(1));
set(AX(1),'YLim',cell2mat(ylimit(1)),'XLim',xlimit,'YColor','k','FontSize',15,...
          'ytick',linspace(tmp(1),tmp(2),4),...
          'xtick',[0,48,96,144,192],...
     'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],'FontSize',15);
    ylabel('%');
    set(Hplot1,'Color','k','LineStyle', '-','LineWidth',2);
 tmp=cell2mat(ylimit(2));   
set(AX(2),'YLim',cell2mat(ylimit(2)),'XLim',xlimit,'YColor','r','FontSize',15,...
           'ytick',linspace(tmp(1),tmp(2),4),...
           'xtick',[0,48,96,144,192],...
     'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],'FontSize',15);
set(get(AX(2),'Ylabel'),'string','%','FontSize',15) 
    set(Hplot2,'Color','r','LineStyle', '-','LineWidth',2);
%set(gca,'xtick',[0,48,96,144,192],...
 %    'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],'FontSize',15);

legend('Elmer','BISICLES');title('(V_{b3}/V_{ASF})x100%');

if 0
 xlimit=[0 204];
 Edata_fwd_asf = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
               cell2mat(area_name(1)),'mb.txt']));
 Edata_control_asf = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_20thsmb_',...
               cell2mat(area_name(1)),'mb.txt']));
 Edata_fwd_b3 = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
               cell2mat(area_name(2)),'mb.txt']));
 Edata_control_b3 = load(cell2mat([E_path_name,'/',base_name_Elmer,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_20thsmb_',...
               cell2mat(area_name(2)),'mb.txt']));    
           
 Bdata_fwd_asf = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
                temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
                cell2mat(area_name(1)),'mb.tab']));     
 Bdata_control_asf = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_20thsmb_',...
               cell2mat(area_name(1)),'mb.tab']));
 Bdata_fwd_b3 = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
                temp_name,'_', cell2mat(beta_name(1)),'_',cell2mat(smb_name(1)),'_',...
                cell2mat(area_name(2)),'mb.tab'])); 
 Bdata_control_b3 = load(cell2mat([B_path_name,'/',base_name_BISICLES,'_',time_name,'_temp',...
               temp_name,'_', cell2mat(beta_name(1)),'_20thsmb_',...
               cell2mat(area_name(2)),'mb.tab']));
           
           
 figure,plot(0:1:monthinyr*(dataed)-1,(Edata_fwd_asf(:,1)-Edata_control_asf(:,1))/1e9,'-k');hold on;
 plot(0:1:monthinyr*(dataed)-1,(Bdata_fwd_asf(:,2)-Bdata_control_asf(:,2))/1e9,'-r');
 xlim(xlimit);
 set(gca,'xtick',[0,48,96,144,192],...
             'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],'FontSize',15);
 ylabel('Volume (km^3)','FontSize',15);
 legend('Elmer','BISICLES');title('ASF: SMB_{elevcorr} - SMB_{20th}')
 
 
 
 figure,plot(0:1:monthinyr*(dataed)-1,(Edata_fwd_b3(:,1)-Edata_control_b3(:,1))/1e9,'-k');hold on;
 plot(0:1:monthinyr*(dataed)-1,(Bdata_fwd_b3(:,2)-Bdata_control_b3(:,2))/1e9,'-r');
 xlim(xlimit);
 set(gca,'xtick',[0,48,96,144,192],...
             'xticklabel',['1995Jan';'1999Jan';'2003Jan';'2007Jan';'2011Jan'],'FontSize',15);
 ylabel('Volume (km^3)','FontSize',15);
 legend('Elmer','BISICLES');title('B3: SMB_{elevcorr} - SMB_{20th}')
end
 