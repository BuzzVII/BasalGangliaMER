using Plots

function h(t, rate, dt, j)
  k = floor(Int,t/dt)
  i = 1
  if k > j[1]
    i = j[1]
  elseif k >= 1
    i = k
  end
  return rate[i]
end

ξ = 1.0
N = 10000
T = 1.0
ve(t) = 14.0 - 4.0*sin(2*30.0*pi*t)
delays = true

#simulation paramters
fs = 24000.0
dt = 1/fs
STN_max_v = 500.0

include("basal_ganglia.jl")
bg_model, ζ = basal_ganglia(ve, delays, ξ)

#Neuron parameters
It = squeeze(readcsv("apcurrent24k.dat"),1)
It = It / maximum(abs(It))
dist = (rand(1,N).^0.3 * 0.001).^-2
scale = dist
scale = scale/maximum(scale)

#the simulation
neuron_superposition = zeros(Int(fs*T))
rates = zeros(5, Int(fs*T))
U = [9.293; -18.4203; 3.857; -4.0643; 16.2714; -194.239; 11.39; -135.967; -1.544; -18.446]
rates[:,1] = [ζ(U[1], θ_p1, k_p1, S_p1) ζ(U[3], θ_p2, k_p2, S_p2) ζ(U[5], θ_d1,  k_d1,  S_d1) ζ(U[7], θ_d2,  k_d2,  S_d2) ζ(U[9], θ_Ϛ, k_Ϛ, S_Ϛ)]'
k1 = zeros(10,1)
k2 = zeros(10,1)
k3 = zeros(10,1)
k4 = zeros(10,1)
refactory = zeros(N,1)
neurons = zeros(N,1)
j = [1]
g = [t -> h(t, rates[ceil(Int, k/2),:], dt, j) for k in 1:10]

for i = 2:floor(Int,fs*T)
    j[1] = i
    #calculate BG rates using RK4
    bg_model(i*dt, U, k1, g)
    bg_model(i*dt + 0.5*dt, U + 0.5*k1*dt, k2, g)
    bg_model(i*dt + 0.5*dt, U + 0.5*k2*dt, k3, g)
    bg_model(i*dt + dt, U + k3*dt, k4, g)
    U += (k1 + 2*k2 + 2*k3 + k4)*dt/6
    rates[:,i] = [ζ(U[1], θ_p1, k_p1, S_p1) ζ(U[3], θ_p2, k_p2, S_p2) ζ(U[5], θ_d1,  k_d1,  S_d1) ζ(U[7], θ_d2,  k_d2,  S_d2) ζ(U[9], θ_Ϛ, k_Ϛ, S_Ϛ)]
    #Generate MER
    P  = dt * rates[5,i]
    refactory = refactory - STN_max_v/fs
    refactory[refactory .< 0] = 0
    neurons = rand(N,1) .<= P
    neurons[refactory .> 0] = 0
    refactory = refactory + neurons
    neuron_superposition[i] = (scale * neurons)[1]
end
MER = conv(neuron_superposition, It)
MER = MER + 0.001*randn(size(MER))
MER = MER - mean(MER)
