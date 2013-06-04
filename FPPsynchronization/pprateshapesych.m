function [t Vout PP comp_time]=pprateshapesych(c,N,groups,neuron_synch_number)
%function [t Vout]=pprateshapesych simulates a patient microelectrode recording
%using a filtered point process with spikes emergant from synchronisation. 
%This is simulation is based on pprateshape, but extends it by adding in
%the ability to have groups of neurons that also fire synchronized based on
%a interval pdf.

tic
start_cpu_time=cputime;
write_file = 1;

h=waitbar(0,'Simulation Progress');

%% Initialize parameters and check arguments passed
if (nargin == 0)
    c=1;
    N=10000;                    
    groups=1;
    neuron_synch_number=100;
end

if (nargin == 1 || N<1)
    N=10000;       
    groups=1;
    neuron_synch_number=100;
end

if (nargin == 2 || groups<1)
    groups=1;
    neuron_synch_number=100;
end

if (nargin == 3 || neuron_synch_number<0)
    neuron_synch_number=100;
end

rate=10;                        %spike rate
lambda=1/rate*gamma(1+1/c);     %calculate scale parameter based on rate and shape
tmax=2;                         %simulation time length
dt=1/24000;                     %time step size
t=0:dt:tmax;                    %create simulation time vector

synchrate=10;                   %synchrony rate
synch_times=zeros(groups,synchrate);
synch_group=zeros(groups,neuron_synch_number);

local_groups=[800 1500 2500 4000 7000 10000];
start_groups=[0 500 1000 2000 3000 5000];
for i=1:groups
    synch_times(i,:)=(tmax*rand(1,synchrate));      %randomly generates synch firing times
    synch_group(i,:)=local_groups(i).*round(rand(1,neuron_synch_number))+start_groups(i);%round(((i-1)+rand(1,neuron_synch_number))*N/6);  
end

%% Create AP waveform
It=dlmread('apcurrent24k.dat');        %Read in current waveform of action potential sin(4000*2*pi*(0:1/24000:30/24000)).*exp(-24000*(0:1/24000:30/24000));%
It=-It./min(It).*250e-9;               %normalize
length_curr=length(It)-1;
epsilon=8.85e-12;                      %Permitivity of free space
sigma=  4.8;                                %conductivity of white matter
rho=10^4*10^6;                          %density of neurons in STN m^-3
r=(3/4*N*rand(N,1)/(pi*rho)).^(1/3);   %create a power law distribution of neuron radii
rsort=sort(r);


weight=ones(size(r));
R3=0.96e3;
C3=2.22e-6;
C2=9.38e-9;
C3=1.56e-6;
C2=9.38e-9;
R4=100e6;
R2N=1./(4*pi*sigma*rsort(length(rsort):-1:1));
R1=2100;
t_impulse=0:1/24000:100/24000;

%prepare voltage and PP storage vector
Vt=zeros(length(t),1).';
Z=Vt;
PP=1;%zeros(N,length(t));

fprintf('simulation initialization complete\n')


% filt=0;

%% Simulate each neuron as a filtered point process
for neuron=1:N
    waitbar(neuron/N,h)
    %initialize vector sizes
    ppwave=zeros(1,length(t));
    pp=ppwave;
    tk=zeros(1,ceil(rate.*max(t)));
    
    %randomly create isi distribution
    isiNonShifted=randraw('weibull',[0,c,lambda],[1,5.*rate.*ceil(max(t))]);
    shift_amount=round(rand*length(isiNonShifted));
    isiShifted=circshift(isiNonShifted,shift_amount);
    shift_amount=round(rand*length(isiShifted));
    isi=circshift(isiShifted,shift_amount);
    
    
    %randomly start the first neuron firing time
    tkindex=2;
    isiindex=2;
    tk(1)=isi(1);
    
    %find then absolute time of each spike time
    if max(max(synch_group == neuron))           %neuron is in synchronized group
        for i=1:groups
%             if i==1
%                 filt=1;
%             end
            if max(synch_group(i,:)==neuron)
                spike_times_temp=sort(synch_times(i,:));
                synch_number=1;
                while tk(tkindex-1)+length_curr*dt<tmax;
                    if synch_number<length(spike_times_temp) && spike_times_temp(synch_number)>tk(tkindex-1) && spike_times_temp(synch_number)<tk(tkindex-1)+isi(isiindex)
                        tk(tkindex)=spike_times_temp(synch_number);   
                        synch_number=synch_number+1;
                        tkindex=tkindex+1;
                    elseif synch_number<length(spike_times_temp) && spike_times_temp(synch_number)==tk(tkindex-1)+isi(isiindex)
                        tk(tkindex)=tk(tkindex-1)+isi(isiindex);
                        tkindex=tkindex+1;
                        isiindex=isiindex+1;
                        synch_number=synch_number+1;
                    else
                        tk(tkindex)=tk(tkindex-1)+isi(isiindex);
                        tkindex=tkindex+1;
                        isiindex=isiindex+1;
                    end
                end
            end
        end
    %find absolute time if no coincidence firing
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
%         if filt
%             [b,a]=butter(8,0.01,'low');
%             If=filter(b,a,It);
%             ppwave(wave_start:wave_end)=ppwave(wave_start:wave_end)+If;   
%         else
            ppwave(wave_start:wave_end)=ppwave(wave_start:wave_end)+It;  
%         end
    end
   

    
   
    %PP(neuron,:)=pp;
    R2=R2N(neuron);
    
    extracellular_impulse_response=-(R4*exp(-(t_impulse*(C2*R1*R2 + C2*R1*R3 + C2*R1*R4 - C3*R1*R3 + C3*R2*R3 + C3*R3*R4))/(2*C2*C3*R1*R3*(R2 + R4))).*(cosh((t_impulse*(C2^2*R1^2*R2^2 + 2*C2^2*R1^2*R2*R3 + 2*C2^2*R1^2*R2*R4 + C2^2*R1^2*R3^2 + 2*C2^2*R1^2*R3*R4 + C2^2*R1^2*R4^2 + 2*C2*C3*R1^2*R2*R3 - 2*C2*C3*R1^2*R3^2 + 2*C2*C3*R1^2*R3*R4 - 2*C2*C3*R1*R2^2*R3 - 2*C2*C3*R1*R2*R3^2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3^2*R4 - 2*C2*C3*R1*R3*R4^2 + C3^2*R1^2*R3^2 - 2*C3^2*R1*R2*R3^2 - 2*C3^2*R1*R3^2*R4 + C3^2*R2^2*R3^2 + 2*C3^2*R2*R3^2*R4 + C3^2*R3^2*R4^2)^(1/2))/(2*C2*C3*R1*R3*(R2 + R4))) + (sinh((t_impulse*(C2^2*R1^2*R2^2 + 2*C2^2*R1^2*R2*R3 + 2*C2^2*R1^2*R2*R4 + C2^2*R1^2*R3^2 + 2*C2^2*R1^2*R3*R4 + C2^2*R1^2*R4^2 + 2*C2*C3*R1^2*R2*R3 - 2*C2*C3*R1^2*R3^2 + 2*C2*C3*R1^2*R3*R4 - 2*C2*C3*R1*R2^2*R3 - 2*C2*C3*R1*R2*R3^2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3^2*R4 - 2*C2*C3*R1*R3*R4^2 + C3^2*R1^2*R3^2 - 2*C3^2*R1*R2*R3^2 - 2*C3^2*R1*R3^2*R4 + C3^2*R2^2*R3^2 + 2*C3^2*R2*R3^2*R4 + C3^2*R3^2*R4^2)^(1/2))/(2*C2*C3*R1*R3*(R2 + R4)))*(C2*R1*R2 - C2*R1*R3 + C2*R1*R4 + C3*R1*R3 - C3*R2*R3 - C3*R3*R4))/(C2^2*R1^2*R2^2 + 2*C2^2*R1^2*R2*R3 + 2*C2^2*R1^2*R2*R4 + C2^2*R1^2*R3^2 + 2*C2^2*R1^2*R3*R4 + C2^2*R1^2*R4^2 + 2*C2*C3*R1^2*R2*R3 - 2*C2*C3*R1^2*R3^2 + 2*C2*C3*R1^2*R3*R4 - 2*C2*C3*R1*R2^2*R3 - 2*C2*C3*R1*R2*R3^2 - 4*C2*C3*R1*R2*R3*R4 - 2*C2*C3*R1*R3^2*R4 - 2*C2*C3*R1*R3*R4^2 + C3^2*R1^2*R3^2 - 2*C3^2*R1*R2*R3^2 - 2*C3^2*R1*R3^2*R4 + C3^2*R2^2*R3^2 + 2*C3^2*R2*R3^2*R4 + C3^2*R3^2*R4^2)^(1/2)))/(C2*(R2 + R4));
        
    electrode_ppwave=conv(ppwave,extracellular_impulse_response,'same');
    

    %convolve delta spikes with action potential to create a voltage spike
    %train, then add to previous neuron voltage spike trains
    %Vw=Vw+weight(neuron).*dt.*fft([It, zeros(1,length(pp)-length(It))]).*fft(pp);
    Vt=Vt+electrode_ppwave;
    
    if ~mod(neuron,1000)
        fprintf('neurons calculated: %i\n',neuron);
        fprintf('time elapsed: %d\n',cputime-start_cpu_time);
    end
    
    
%     filt=0;
  
    
end

close(h)
fprintf('neuron contribution to MER complete\n'); %neuron contribution to MER complete

%% MER equipment effects

%change the MER probe voltage to a time series (faster for complete voltage
%signal than doing for each neuron individually).
%Vt=(ifft(Vw))/length(Vw);

%mean zero the signal
Vt=Vt-mean(Vt);

%set high and low pass filter corners
flow=5000;
fhigh=500;

%apply phase preserving filters to signal
% [b,a]=butter(9,flow/24000,'low');
% Vf=filter(b,a,Vt);
% [b,a]=butter(2,fhigh/24000,'high');
% Vf=filter(b,a,Vf);

%set electrode parameters to calculate random noise
R=1e6;
kB=1.38e-23;
T=37+273;
Nstd=4.*kB.*T.*R.*0;
Ntherm=Nstd.*randn(size(t));

%add random noise
Vout=Vt+Ntherm;

% if write_file
%     writeMERsimulation(N,'weibull',0,c,lambda,Vout,PP,'whitenoise',num2str(Nstd),1,N);
%     fprintf('simulation written to file ppsim%s.dat\nTotal run time: %d\n',1,cputime-start_cpu_time)
% end

comp_time=toc;

fprintf('simulation complete in %d seconds\n',comp_time) %simulation complete


