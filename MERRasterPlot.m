function MERRasterPlot(timevec,DataMER,DataPP,N)%ppfile,fix_stim)


 %[header,DataMER,ppheader,DataPP]=readMERsimulation(ppfile);
 
neuron_properties=struct('colour',[]);
% 
% if fix_stim<=0
%     load('StimulusData_th0.mat');
%     patient_time=(0:(length(stim_struct.x)-1))*1000/stim_struct.fs;
%     patient_signal=stim_struct.x/10;
% else
%     load('FixationData_th0.mat');
%     patient_time=(0:(length(fix_struct.x)-1))*1000/fix_struct.fs;
%     patient_signal=fix_struct.x/10;
% end

for marker_colour = 1:N
		      neuron_properties(marker_colour).colour=[1-marker_colour/N 0 1-(N-marker_colour)/N];
end

 figure(1)
 clf
 subplot(2,1,1)
 hold on
 
 for timestep=1:length(timevec)%1*ppheader.fs+1:4.1*ppheader.fs
    time=timevec(timestep);%(timestep-1)/ppheader.fs;
%     neuron_number = 1;
%     neurons_found = 0;
%     N=[];
%     coincidence=0;
%     coincidence_count=0;
%      for power=0:N-1
%         if mod(DataPP(timestep)-neurons_found,2^power)
%             N(neuron_number)=power;
%             neuron_number=neuron_number+1;
%             coincidence=coincidence+1;
%             neurons_found=neurons_found+2^(power-1);
%         end
%      end
%      if coincidence > 3
%          coincidence_count=coincidence_count+1;
%          coincidences(coincidence_count)=time;
%      end
     for neuron=1:N-1
		  %plot(time,N(neuron),'LineWidth',1,'LineStyle','--','color',[1 1 1],'Marker','square','MarkerFaceColor',neuron_properties(N(neuron)).colour,'MarkerEdgeColor',neuron_properties(N(neuron)).colour,'MarkerSize',3)
          if  DataPP(neuron,timestep)
            plot(time,neuron,'LineWidth',1,'LineStyle','--','color',[1 1 1],'Marker','square','MarkerFaceColor',neuron_properties(neuron).colour,'MarkerEdgeColor',neuron_properties(neuron).colour,'MarkerSize',0.05)
            %pause(0.1);
          end
     end
 end
 
 %ylim([0 65]) 
 %xlim([1 2])
 set(gca,'box','on','fontsize',16)%,'OuterPosition',[0 0.4 1 0.5])%,'XTickLabelMode','manual','XTickLabel',[],'YTickMode','manual','YTick',[2; 4; 6; 8])
 ylabel('neuron','fontsize',16)
 
 
 subplot(2,1,2)
 plot(timevec,DataMER); %(1*header.fs:2*header.fs)/header.fs*1000, -DataMER((1*header.fs):(2*header.fs)));
 %xlim([1000 2000])


 set(gca,'box','on','fontsize',16)%,'OuterPosition',[0 0 1 0.5])
 ylabel('Voltage(mV)','fontsize',16)%,'Position',[3.99 0 17.3])
 
%  subplot(3,1,3)
%  plot(patient_time,patient_signal,'g');
% xlim([0 1000])
% set(gca,'box','on','fontsize',16) 
%  xlabel('time (ms)','fontsize',16)
% ylabel('Voltage(mV)','fontsize',16)

 set(gcf,'color',[1 1 1])
