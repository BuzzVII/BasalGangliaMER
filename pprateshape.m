function [t Vout comp_time]=pprateshape(c,N,synchronization,fix_stim)
%function [t Vout]=pprateshape simulates a patient microelectrode recording
%using a filtered point process. This is done by simulating N neuron spike 
%times and then applying extracellular filtering to create the effective
%voltage the electrode 'sees' for each neuron. This signal then undergoes
%processing that matches the equipment effects from a real recording 
%(including electrode transfer function, hardware filtering, electronic
%noise). The simulation is then written to a pp simulation data file.

tic
start_cpu_time=cputime;

%% Initialize parameters and check arguments passed
if (nargin == 1 || N<1)
    N=10000;                    %Total number of neurons
end

if nargin == 2
    synchronization=0;          %determines if any neurons have coupled firing times
    fix_stim=0;
elseif nargin == 3
    fix_stim=0;
end

if synchronization
    neuron_synch_number=100;
    if fix_stim
        load('FixationData_th0.mat');
        thresholdTime=(1:length(fix_struct.thresholdSignal)-1)/fix_struct.fs;
        thresholdSpikes=thresholdTime(fix_struct.thresholdSignal>1);
        synch_times=thresholdSpikes+1;%(4+rand(1,100));      %randomly generates 100 synch firing times
        synch_chance=boolean(round(19/32*rand(neuron_synch_number,length(synch_times))));  %chance of each neuron firing at synch times   
    else
        load('StimulusData_th0.mat');
        thresholdTime=(1:length(stim_struct.thresholdSignal)-1)/stim_struct.fs;
        thresholdSpikes=thresholdTime(stim_struct.thresholdSignal>1);
        synch_times=thresholdSpikes+1;%(4+rand(1,100));      %randomly generates 100 synch firing times
        synch_chance=boolean(round(5/8*rand(neuron_synch_number,length(synch_times))));  %chance of each neuron firing at synch times   
    end
end

rate=30;                        %spike rate
%c=100;                         %weibull shape parameter
lambda=1/rate*gamma(1+1/c);     %calculate scale parameter based on rate and shape
tmax=20;                        %simulation time length
dt=1/24000;                     %time step size
t=0:dt:tmax;                    %create simulation time vector


%% Create AP waveform
It=dlmread('apcurrent24k.dat');        %Read in current waveform of action potential
It=-It./min(It).*250e-9;               %normalize
length_curr=length(It)-1;
epsilon=8.85e-12;                      %Permitivity of free space
rho=10^5*10^6;         	   	  	%density of neurons in STN m^-3
r=(3/4*N*rand(N,1)/(pi*rho)).^(1/3);   %create a power law distribution of neuron radii
rsort=sort(r);

%weight=1./(4*pi*epsilon*rsort);          %current to voltage coefficient
weight=ones(size(r));
R3=0.96e3;
C3=2.22e-6;
C2=9.38e-9;
C3=1.56e-6;
C2=9.38e-9;
R4=100e6;
R2N=1./(4*pi*epsilon*rsort(length(rsort):-1:1));
R1=2100;
t_impulse=0:1/24000:100/24000;


% figure(22)
% hist(r);
%Zw=dlmread('Zwinterp.dat');            %Extracellular filter function
%

%prepare voltage and PP storage vector
Vt=zeros(length(t),1).';
Z=Vt;
PP8=Z;





fprintf('simulation initialization complete\n')    %initialization complete


%% Simulate each neuron as a filtered point process
for neuron=1:N
    
    %initialize vector sizes
    ppwave=zeros(1,length(t));
    pp=ppwave;
    tk=zeros(1,ceil(rate.*max(t)));
    
    %randomly create isi distribution
    isiNonShifted=randraw('weibull',[0,c,lambda],[1,5.*rate.*round(max(t))]);
    shift_amount=round(rand*length(isiNonShifted));
    isiShifted=circshift(isiNonShifted,shift_amount);
    shift_amount=round(rand*length(isiShifted));
    isi=circshift(isiShifted,shift_amount);
    
    
    %randomly start the first neuron firing time
    tkindex=2;
    isiindex=2;
    tk(1)=isi(1);%+1/rate*rand(1,1);%randraw('exp',rate,1);
    
    %find then absolute time of each spike time
    %add coincidence firing times for the first 64 neurons     
    if synchronization && neuron<=neuron_synch_number 
        spike_times_temp=sort(synch_times(synch_chance(neuron,:)));
        synch_number=1;
        %keyboard;
        while tk(tkindex-1)+length_curr*dt<tmax;
            if synch_number<length(spike_times_temp) && spike_times_temp(synch_number)>tk(tkindex-1) && spike_times_temp(synch_number)<tk(tkindex-1)+isi(isiindex)
                tk(tkindex)=spike_times_temp(synch_number);   
                synch_number=synch_number+1;
                tkindex=tkindex+1;
                %fprintf('synch time written for neuron %i at time:%d\n',neuron,spike_times_temp(synch_number));
            elseif synch_number<length(spike_times_temp) && spike_times_temp(synch_number)==tk(tkindex-1)+isi(isiindex)
                tk(tkindex)=tk(tkindex-1)+isi(isiindex);
                tkindex=tkindex+1;
                isiindex=isiindex+1;
                synch_number=synch_number+1;
                %fprintf('synch time already exists for neuron:%i\n',neuron);
            else
                tk(tkindex)=tk(tkindex-1)+isi(isiindex);
                tkindex=tkindex+1;
                isiindex=isiindex+1;
            end
        end
    %find absolut time if no coincidence firing
    else
        while tk(tkindex-1)+length_curr*dt<tmax;
            tk(tkindex)=tk(tkindex-1)+isi(isiindex);
            tkindex=tkindex+1;
            isiindex=isiindex+1;
        end
    end
  
    
    
    
    %set a delta spike at each spike timing
    for i=1:tkindex-2
        wave_start=round(tk(i)/dt)+1;
        wave_end=wave_start+length_curr;
        pp(wave_start)=1;
        ppwave(wave_start:wave_end)=ppwave(wave_start:wave_end)+It;
    end
    
    
%     shift_amount=round(rand*length(ppwave));
%     ppwave=circshift(ppwave,shift_amount);
%     pp=circshift(ppwave,shift_amount);
    
    R2=R2N(neuron);
    
    extracellular_impulse_response=-(R4*exp(-(t_impulse*(C2*R1*R2 + C2*R1*R3 + C2*R1*R4 - C3*R1*R3 + C3*R2*R3 + C3*R3*R4))/(2*C2*C3*R1*R3*(R2 + R4))).*(cosh((t_impulse*(C2^2*R1^2*R2^2 + 2*C2^2*R1^2*R2*R3 + 2*C2^2*R1^2*R2*R4 + C2^2*R1^2*R3^2 + 2*C2^2*R1^2*R3*R4 + C2^2*R1^2*R4^2 + 2*C2*C3*R1^2*R2*R3 - 2*C2*C3*R1^2*R3^2 + 2*C2*C3*R1^2*R3*R4 - 2*C2*C3*R1*R2^2*R3 - 2*C2*C3*R1*R2*R3^2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3^2*R4 - 2*C2*C3*R1*R3*R4^2 + C3^2*R1^2*R3^2 - 2*C3^2*R1*R2*R3^2 - 2*C3^2*R1*R3^2*R4 + C3^2*R2^2*R3^2 + 2*C3^2*R2*R3^2*R4 + C3^2*R3^2*R4^2)^(1/2))/(2*C2*C3*R1*R3*(R2 + R4))) + (sinh((t_impulse*(C2^2*R1^2*R2^2 + 2*C2^2*R1^2*R2*R3 + 2*C2^2*R1^2*R2*R4 + C2^2*R1^2*R3^2 + 2*C2^2*R1^2*R3*R4 + C2^2*R1^2*R4^2 + 2*C2*C3*R1^2*R2*R3 - 2*C2*C3*R1^2*R3^2 + 2*C2*C3*R1^2*R3*R4 - 2*C2*C3*R1*R2^2*R3 - 2*C2*C3*R1*R2*R3^2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3^2*R4 - 2*C2*C3*R1*R3*R4^2 + C3^2*R1^2*R3^2 - 2*C3^2*R1*R2*R3^2 - 2*C3^2*R1*R3^2*R4 + C3^2*R2^2*R3^2 + 2*C3^2*R2*R3^2*R4 + C3^2*R3^2*R4^2)^(1/2))/(2*C2*C3*R1*R3*(R2 + R4)))*(C2*R1*R2 - C2*R1*R3 + C2*R1*R4 + C3*R1*R3 - C3*R2*R3 - C3*R3*R4))/(C2^2*R1^2*R2^2 + 2*C2^2*R1^2*R2*R3 + 2*C2^2*R1^2*R2*R4 + C2^2*R1^2*R3^2 + 2*C2^2*R1^2*R3*R4 + C2^2*R1^2*R4^2 + 2*C2*C3*R1^2*R2*R3 - 2*C2*C3*R1^2*R3^2 + 2*C2*C3*R1^2*R3*R4 - 2*C2*C3*R1*R2^2*R3 - 2*C2*C3*R1*R2*R3^2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3^2*R4 - 2*C2*C3*R1*R3*R4^2 + C3^2*R1^2*R3^2 - 2*C3^2*R1*R2*R3^2 - 2*C3^2*R1*R3^2*R4 + C3^2*R2^2*R3^2 + 2*C3^2*R2*R3^2*R4 + C3^2*R3^2*R4^2)^(1/2)))/(C2*(R2 + R4));
    electrode_ppwave=conv(ppwave,extracellular_impulse_response,'same');
    
    if neuron<64
        PP8=PP8+2^(63-neuron)*pp;
    end

    %convolve delta spikes with action potential to create a voltage spike
    %train, then add to previous neuron voltage spike trains
    %Vw=Vw+weight(neuron).*dt.*fft([It, zeros(1,length(pp)-length(It))]).*fft(pp);
    Vt=Vt+electrode_ppwave;
    
    if ~mod(neuron,200)
        fprintf('neurons calculated: %i\n',neuron);
        fprintf('time elapsed: %d\n',cputime-start_cpu_time);
        
    end
end

fprintf('neuron contribution to MER complete\n'); %neuron contribution to MER complete

%% MER equipment effects

%change the MER probe voltage to a time series (faster for complete voltage
%signal than doing for each neuron individually).
%Vt=(ifft(Vw))/length(Vw);

%mean zero the signal
Vt=Vt-mean(Vt);

%set high and low pass filter corners
flow=10000;
fhigh=500;

%apply phase preserving filters to signal
[b,a]=butter(9,flow/24000,'low');
Vf=filter(b,a,Vt);
% [b,a]=butter(0,fhigh/24000,'high');
% Vf=filter(b,a,Vf);

%set electrode parameters to calculate random noise
R=1e6;
kB=1.38e-23;
T=37+273;
Nstd=4.*kB.*T.*R.*0;
Ntherm=Nstd.*randn(size(t));

%add random noise
Vout=Vf+Ntherm;

% figure(2)
% plot(t,pp);
% xlabel('t (s)');title('Single neuron spike train');

fprintf('simulation complete\n') %simulation complete


 Nid=length(dir('h:\Data\ppsim*'))+1;                         %Output file ID
% %% Write simulation to data file
 writeMERsimulation(N,'weibull',0,c,lambda,Vout,PP8,'whitenoise',num2str(Nstd),Nid,N);
% 
 fprintf('simulation written to file ppsim%s.dat\nTotal run time: %d\n',num2str(Nid),cputime-start_cpu_time) %simulation written to file
% 
%% Plots
%figure(1)
%plot(t,Vout);
%xlabel('t (s)');ylabel('V (V)');title('MER simulation');
%print -depsc MERvoltage

comp_time=toc
