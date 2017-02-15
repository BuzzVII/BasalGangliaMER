using Plots

function observability(C,A)
  D = C
  for n  = 1:(size(A)[1] - 1)
    D = [D; C*A^n]
  end
  rank(D)
end

ξ = 1.0
N = 10000
T = 1.0
ve(t) = 14.0 - 4.0*sin(2*30.0*pi*t) + 4*rand()
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
dU = zeros(10,1)
refactory = zeros(N,1)
neurons = zeros(N,1)
C = zeros(10,10)
C[9,9] = 1; C[10,10] = 1
obs = zeros(1, Int(fs*T))
A = [[0        1         0              0                  0        0       0            0           0           0]
     [(α*β) (β + α) (α*β)*ν[2,1] (α*β)*ν[2,1]     (α*β)*ν[3,1] (α*β)*ν[3,1] 0            0     (α*β)*ν[5,1] (α*β)*ν[5,1]]
     [0        0         0              1                  0        0       0            0           0           0]
     [0        0 (α*β)+(α*β)*ν[2,2] (β + α)+(α*β)*ν[2,2]   0        0  (α*β)*ν[4,2] (α*β)*ν[4,2] (α*β)*ν[5,2] (α*β)*ν[5,2]]
     [0        0         0              0                  0        1       0            0           0           0]
     [0        0         0              0                (α*β)   (β + α)    0            0           0           0]
     [0        0         0              0                  0        0       0            1           0           0]
     [0        0         0              0                  0        0     (α*β)        (β + α)       0           0]
     [0        0         0              0                  0        0       0            0           0           1]
     [0        0    (α*β)*ν[2,5]    (α*β)*ν[2,5]           0        0       0            0         (α*β)     (β + α)]]

for i = 2:Int(fs*T)
    #calculate BG rates
    bg_model(i*dt,U,dU)
    U += dU*dt
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
