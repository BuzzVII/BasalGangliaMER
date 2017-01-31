function fpp(DA, N, T)
  #             GPi    GPe  D1   D2 STN  FS  Ctx         Tha  TRN
  weights =  [[   0   -0.3  -1    0 0.9   0    0          0    0]
             [    0      0   0   -1 0.9   0    0          0    0]
             [    0 (1+DA)   0    0   0  -1 (1+DA)*0.5    0    0]
             [    0 (1-DA)   0    0   0  -1 (1-DA)*0.5    0    0]
             [    0     -1   0    0   0   0    0.5        0    0]
             [    0      1   0    0   0   0    0          0    0]
             [    0      0   0    0   0   0    0        0.6    0]
             [ 0.18      0   0    0   0   0    0.6        0 0.35]
             [    0      0   0    0   0   0    1       0.35    0]];
  #in ms
  delays =  [[    0      0   0    5   2   0    0          0    0]
             [    0      0   0    5   2   0    0          0    0]
             [    0      0   0    0   0   0    10         0    0]
             [    0      0   0    0   0   0    10         0    0]
             [    4      4   0    0   0   0    3          0    0]
             [    0      0   0    0   0   0    0          0    0]
             [    0      0   0    0   0   0    0          0    0]
             [    0      0   0    0   0   0    0          0    0]
             [    0      0   0    0   0   0    0          0    0]];
  # time constant for membrane potentials
  tau = [0.02, 0.02, 0.02, 0.02, 0.005, 0.005, 0.08, 0.005, 0.005];
  excitatory = [-1, -1, -1, -1, 1, 1, 1, 1, 1];
  current = zeros(9,9);
  DBS = zeros(9,1);
  R = 100 * ones(9,1);

  #sigmoid parameters
  max_rate = 0.5*ones(9,1);
  slope = ones(9,1);
  max_potential = zeros(9,1);

  #simulation paramters
  fs = 24000;
  dt = 1/fs;
  STN_max_v = 111;

  #Neuron parameters
  It = squeeze(readcsv("apcurrent24k.dat"),1);
  It = It / maximum(abs(It));
  dist = (rand(1,N).^0.3 * 0.001).^-2;
  scale = dist;
  scale = scale/maximum(scale);

  #the simulation
  neuron_superposition = zeros(Int(fs*T));
  rates = zeros(9, Int(fs*T));
  rates[7,1] = 1;
  refactory = zeros(N,1);
  delays = delays * Int(fs / 1000);
  neurons = zeros(N,1);
  for i = 2:Int(fs*T)
      #calculate BG rates
      current_delay = max(min((i - 1) - delays, i-1)-1, 0) * 9 + repmat((1:9)',9,1);
      rate = rates[current_delay];
      #solve the new synaptic currents using: tau*dJ/dt=-J+v
      current = current + dt./tau' .* (rate - current);
      #sum all the synaptic currents into the neuron weighted by a gain term G: I=sum(G*j)
      synaptic_current = weights * current';
      synaptic_current = diag(synaptic_current) + DBS;
      u = R .* synaptic_current + excitatory .* rates[:,i-1];
      sigma_rate = 2 * max_rate ./ (1 + exp(slope .* (max_potential - u)));
      rand_input = zeros(9,1);
      rand_input[7] = 1;
      #solve the firing rate for the structure using the synaptic inputs: tau*dv/dt=-v+S[R*I+/-v]
      rates[:,i] = rates[:,i-1] + dt ./ tau .* ( sigma_rate - rates[:,i-1] + rand_input );
      rates[rates[:,i] .> 1, i] = 1;
      rates[rates[:,i] .< 0, i] = 0;

      #Generate MER
      P  = dt * STN_max_v * rates[5,i];
      refactory = refactory - STN_max_v/fs;
      refactory[refactory .< 0] = 0;
      neurons = rand(N,1) .<= P;
      neurons[refactory .> 0] = 0;
      refactory = refactory + neurons;
      neuron_superposition[i] = (scale * neurons)[1];
  end
  MER = conv(neuron_superposition, It);
  MER = MER + 0.001*randn(size(MER));
  MER = MER - mean(MER);
  return MER, rates
end
