function generate_timing(T, fs, sol, N)
  dt = 1/fs
  STN_max_ν = 500.0
  dist = (rand(1,N).^0.3 * 0.001).^-2
  scale = dist
  scale = scale/maximum(scale)
  N_max = floor(Int,fs*T)
  neuron_superposition = zeros(N_max)
  refactory = zeros(N,1)
  neurons = zeros(N,1)
  for i = 0:(N_max - 1)
      P  = dt * ζ(sol(i*dt)[9], θ_Ϛ, k_Ϛ, S_Ϛ)
      refactory = refactory - STN_max_ν/fs
      refactory[refactory .< 0] = 0
      neurons = rand(N,1) .<= P
      neurons[refactory .> 0] = 0
      refactory = refactory + neurons
      neuron_superposition[i+1] = (scale * neurons)[1]
  end
  return neuron_superposition
end
