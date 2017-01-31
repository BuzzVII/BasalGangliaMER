rng(22011985)
%Neural mass model parameters
DA = 0;
%             GPi    GPe  D1   D2 STN  FS  Ctx         Tha  TRN  
weights =	[   0,  -0.3, -1,   0,0.9,  0,   0	    ,   0,   0;
                0,     0,  0,  -1,0.9,  0,   0	    ,   0,   0;
                0,(1+DA),  0,   0,  0, -1,(1+DA)*0.5,   0,   0;
                0,(1-DA),  0,   0,  0, -1,(1-DA)*0.5,   0,   0;
                0,    -1,  0,   0,  0,  0,   0.5    ,   0,   0;
                0,     1,  0,   0,  0,  0,   0	    ,   0,   0;
                0,     0,  0,   0,  0,  0,   0 	    , 0.6,   0;
             0.18,     0,  0,   0,  0,  0,   0.6    ,   0,0.35;
                0,     0,  0,   0,  0,  0,   1	    ,0.35,   0];
%in ms
delays =    [   0,     0,  0,   5,  2,  0,   0	  ,   0,   0;
                0,     0,  0,   5,  2,  0,   0	  ,   0,   0;
                0,     0,  0,   0,  0,  0,   10   ,   0,   0;
                0,     0,  0,   0,  0,  0,   10   ,   0,   0;
                4,     4,  0,   0,  0,  0,   3 	  ,   0,   0;
                0,     0,  0,   0,  0,  0,   0	  ,   0,   0;
                0,     0,  0,   0,  0,  0,   0 	  ,   0,   0;
                0,     0,  0,   0,  0,  0,   0	  ,   0,   0;
                0,     0,  0,   0,  0,  0,   0	  ,   0,   0];
% time constant for membrane potentials
tau = [0.02, 0.02, 0.02, 0.02, 0.005, 0.005, 0.08, 0.005, 0.005]';
excitatory = [-1, -1, -1, -1, 1, 1, 1, 1, 1]';
current = zeros(9);
DBS = zeros(9,1);
R = 100 * ones(9,1);

%sigmoid parameters
max_rate = 0.5*ones(9,1);
slope = ones(9,1);
max_potential = zeros(9,1);

%simulation paramters
fs = 24000;
dt = 1/fs;
v  = 30;
STN_max_v = 111;
T  = 1.0;
N  = 1;

%Neuron parameters
It = dlmread('apcurrent24k.dat');
It = It / max(It);
scale = (rand(1,N).^0.3 * 0.001).^-2;
scale = scale/max(scale);

%the simulation
tic
neuron_superposition = zeros(1, fs*T);
rates = rand(9,1);
refactory = zeros(N,1);
delays = delays * fs/1000; %convert into simulation steps
for i = 2:fs*T
    %calculate BG rates
    current_delay = max(min((i - 1) - delays, i-1)-1, 0) * 9 + repmat(1:9,9,1);
    rate = rates(current_delay);
    %solve the new synaptic currents using: tau*dJ/dt=-J+v
    current = current + bsxfun(@times, dt./tau', rate - current);
	%sum all the synaptic currents into the neuron weighted by a gain term G: I=sum(G*j)
    synaptic_current = weights * current';
    synaptic_current = diag(synaptic_current) + DBS;
    u = R .* synaptic_current + excitatory .* rates(:,i-1);
    sigma_rate = 2 * max_rate ./ (1 + exp(slope .* (max_potential - u)));
    rand_input = zeros(9,1);
    rand_input(1:9) = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1]' .* rand(9,1);
    %solve the firing rate for the structure using the synaptic inputs: tau*dv/dt=-v+S[R*I+/-v]
    rates(:,i) = rates(:,i-1) + dt ./ tau .* ( sigma_rate - rates(:,i-1) + rand_input );
	rates(rates(:,i) > 1, i) = 1;
	rates(rates(:,i) < 0, i) = 0;
    
    %Generate MER
    P  = dt * STN_max_v * rates(5,i);
    refactory = refactory - STN_max_v/fs;
    refactory(refactory < 0) = 0;
    neurons = rand(N,1) <= P;
    neurons(refactory > 0) = 0;
    refactory = refactory + neurons;
    neuron_superposition(i) = scale * neurons;
end
MER = conv(neuron_superposition, It, 'valid');
MER = MER + 0.001*randn(size(MER));
MER = MER - mean(MER);
toc

%Plots and shit
figure(1)
plot((0:size(MER,2)-1)*dt, MER)
figure(2)
[P,F] = pwelch(MER,2^14,2^10,2^14,fs);
loglog(F,P)
xlabel('frequency (Hz)')
ylabel('Power')
axis tight
figure(3)
[S,F,Ts] = spectrogram(MER, 2^12, 2^10, 2^12, fs);
surf(Ts,F,abs(S),'EdgeColor','none','LineStyle','none','FaceLighting','phong')
caxis([0, 70])
set(gca, 'ZScale', 'log')
set(gca, 'YScale', 'log')
xlabel('time (s)')
ylabel('frequency (Hz)')
zlabel('Amplitude')
axis tight
figure(4)
hold on;
plot(0:dt:T-dt, rates(5,:))