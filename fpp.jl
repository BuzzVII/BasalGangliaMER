using Plots


# spikes/s
const S_D1 = 300.0; const S_D2 = 300.0; const S_GPe = 400.0; const S_GPi = 300.0; const S_STN = 500.0
# mV^-1
const k_D1 = 0.3; const k_D2 = 0.3; const k_GPe = 0.2;const k_GPi = 0.2; const k_STN = 0.2
# mV
const vh_D1 = 27; const vh_D2 = 27; const vh_GPe = 14; const vh_GPi = 12; const vh_STN = 18.5
#A mV, D ms^-1,
const A_A = 2*e;   const D_A = 1/4.0*e
const A_N = 0.1*e; const D_N = 1/100.0*e
const A_G = 2*e;   const D_G = 1/6.0*e
#              GPi    GPe  D1   D2  STN  Ctx
const C =  [[   0.0    5   22    0   3    0  ]
            [    0     1    0   33   3    0  ]
            [    0     0    0    0   0   50  ]
            [    0     0    0    0   0   50  ]
            [    0    10    0    0   2    2  ]
            [    0     0    0    0   0    0  ]]'
# delays ms
const τ_1 = 1.0; const τ_3 = 3.0; const τ_4 = 4.0
# Cortical firing Hz
const μ = 3.0; const σ2 = 1.0

function sigmoid(v, vh, k, S_max)
  S = S_max/(1+exp(k*(vh - v)))
end

function bg_model(t,u,du)
  du[1] = u[2]
  du[2] = (A_A*D_A + A_N*D_N)*C[5,1]*sigmoid(u[9], vh_STN, k_STN, S_STN) -
          A_G*D_G*C[2,1]*sigmoid(u[3], vh_GPe, k_GPe, S_GPe) -
          A_G*D_G*C[3,1]*sigmoid(u[5], vh_D1,  k_D1,  S_D1) -
          2*(D_A+D_N-2*D_G)*u[2] - (D_A^2+D_N^2-2*D_G^2)*u[1]
  du[3] = u[4]
  du[4] = (A_A*D_A + A_N*D_N)*C[5,2]*sigmoid(u[9], vh_STN, k_STN, S_STN) -
          A_G*D_G*C[4,2]*sigmoid(u[7], vh_D2,  k_D2,  S_D2) -
          A_G*D_G*C[2,2]*sigmoid(u[3], vh_GPe, k_GPe, S_GPe) -
          2*(D_A+D_N-2*D_G)*u[4] - (D_A^2+D_N^2-2*D_G^2)*u[3]
  du[5] = u[6]
  du[6] = (A_A*D_A + A_N*D_N)*C[6,3]*u[11] - 2*(D_A+D_N)*u[6] - (D_A^2+D_N^2)*u[5]
  du[7] = u[8]
  du[8] = (A_A*D_A + A_N*D_N)*C[6,4]*u[11] - 2*(D_A+D_N)*u[8] - (D_A^2+D_N^2)*u[7]
  du[9] = u[10]
  du[10] = (A_A*D_A + A_N*D_N)*C[6,5]*u[11] +
           (A_A*D_A + A_N*D_N)*C[5,5]*sigmoid(u[9], vh_STN, k_STN, S_STN) -
           A_G*D_G*C[2,5]*sigmoid(u[3], vh_GPe, k_GPe, S_GPe) -
           2*(2*D_A+2*D_N-D_G)*u[10] - (2*D_A^2+2*D_N^2-D_G^2)*u[9]
  du[11] = 0;#-u[11] + sqrt(σ2) * randn() + μ; #not sure if this is right for the cortical input
  return du
end

function fpp(DA, N, T)
  #simulation paramters
  fs = 24000
  dt = 1/fs
  STN_max_v = 111

  #Neuron parameters
  It = squeeze(readcsv("apcurrent24k.dat"),1)
  It = It / maximum(abs(It))
  dist = (rand(1,N).^0.3 * 0.001).^-2
  scale = dist
  scale = scale/maximum(scale)

  #the simulation
  neuron_superposition = zeros(Int(fs*T))
  rates = zeros(6, Int(fs*T))
  U = zeros(11,1)
  U[[1,3,5,7,9,11]]=[vh_GPe, vh_GPi, vh_D1, vh_D2, vh_STN, 0]
  rates[:,1] = U[[1,3,5,7,9,11]]
  dU = zeros(11,1)
  refactory = zeros(N,1)
  neurons = zeros(N,1)
  for i = 2:Int(fs*T)
      #calculate BG rates
      dU = bg_model(i*dt,U,dU)
      U += dU*dt
      rates[:,i] = U[[1,3,5,7,9,11]]
      #Generate MER
      P  = dt * STN_max_v * sigmoid(rates[5,i], vh_STN, k_STN, S_STN)
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
  return MER, rates
end

MER, rates = fpp(1.0,1,1.0)
plot(rates[6,:])
